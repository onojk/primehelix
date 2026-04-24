#!/usr/bin/env python3
"""
primehelix_viz.py — Interactive 3D semiprime structure visualizer

Requires:  pip install pyglet        (tested on pyglet 2.x)
Run:       python primehelix_viz.py
           python primehelix_viz.py --limit 500000 --n 20000

Controls:
  Left drag          rotate view
  Mouse wheel        zoom
  Space              pause / resume auto-rotation
  [ / ]              fewer / more semiprimes shown
  - / =              decrease / increase helix radius
  , / .              decrease / increase helix pitch
  9 / 0              slower / faster rotation
  R                  reset view
  S                  save screenshot
  Esc                quit
"""
from __future__ import annotations

import argparse
import ctypes
import math
import sys
import time
from dataclasses import dataclass
from typing import Dict, List, Tuple

try:
    import pyglet
    from pyglet import gl
    from pyglet.graphics.shader import Shader, ShaderProgram
    import pyglet.math as pm
except ImportError:
    sys.exit("pyglet not found — run: pip install pyglet")

# ── Defaults ──────────────────────────────────────────────────────────────────
DEFAULT_LIMIT        = 1_000_000
DEFAULT_N            = 20_000
DEFAULT_RADIUS       = 42.0
DEFAULT_PITCH        = 0.35        # vertical spacing per lopsided index
DEFAULT_ANGLE_STEP   = 0.24        # radians per lopsided index (≈ golden angle / 2π * twist)
DEFAULT_SPEED        = 0.25        # auto-rotation speed (turns per second * factor)
DEFAULT_POINT_SIZE   = 3.5
BG                   = (0.025, 0.028, 0.038, 1.0)

# Residue-family colors (RGBA float 0-1)
FAM_COLORS: Dict[str, Tuple[float, float, float]] = {
    "1x1":      (0.55, 0.20, 0.85),   # purple
    "1x3":      (0.20, 0.82, 0.35),   # green
    "3x3":      (1.00, 0.42, 0.28),   # coral
    "even":     (0.62, 0.62, 0.62),   # gray
    "other":    (0.82, 0.82, 0.82),
    "balanced": (1.00, 0.75, 0.10),   # amber  (axis / non-lopsided)
}

# ── Pure-Python number theory ─────────────────────────────────────────────────

def _is_prime(n: int) -> bool:
    if n < 2: return False
    if n < 4: return True
    if n % 2 == 0 or n % 3 == 0: return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0: return False
        i += 6
    return True


def prime_sieve(limit: int) -> List[int]:
    if limit < 2:
        return []
    flags = bytearray(b"\x01") * (limit + 1)
    flags[0] = flags[1] = 0
    for i in range(2, math.isqrt(limit) + 1):
        if flags[i]:
            flags[i * i :: i] = b"\x00" * len(flags[i * i :: i])
    return [i for i, f in enumerate(flags) if f]


def residue_family(p: int, q: int) -> str:
    if p == 2 or q == 2:
        return "even"
    a, b = p % 4, q % 4
    if a == 1 and b == 1:   return "1x1"
    if a == 3 and b == 3:   return "3x3"
    if {a, b} == {1, 3}:    return "1x3"
    return "other"


def _is_lopsided(p: int, q: int, bit_gap: int = 8) -> bool:
    return (q.bit_length() - 1) - (p.bit_length() - 1) > bit_gap


@dataclass(frozen=True)
class SP:
    n: int
    p: int
    q: int
    lopsided: bool
    family: str


def generate_semiprimes(limit: int, bit_gap: int = 8) -> List[SP]:
    """All semiprimes n=p·q ≤ limit, p≤q both prime, sorted by n."""
    primes = prime_sieve(limit // 2 + 1)
    out: List[SP] = []
    for i, p in enumerate(primes):
        if p * p > limit:
            break
        for q in primes[i:]:
            n = p * q
            if n > limit:
                break
            out.append(SP(n=n, p=p, q=q,
                          lopsided=_is_lopsided(p, q, bit_gap),
                          family=residue_family(p, q)))
    out.sort(key=lambda s: s.n)
    return out

# ── Geometry ──────────────────────────────────────────────────────────────────

def build_geometry(
    sps: List[SP],
    radius: float,
    pitch: float,
    angle_step: float,
) -> Tuple[List[float], List[float], Dict]:
    """
    Returns (helix_verts, axis_verts, stats).
    Each vert list is flat [x,y,z,r,g,b, ...].
    helix_verts: lopsided semiprimes wound on a helix.
    axis_verts:  non-lopsided placed on the y-axis (x=0, z=0).
    """
    helix: List[float] = []
    axis:  List[float] = []
    fam_counts: Dict[str, int] = {}
    lop_i = bal_i = 0

    n_lopsided = sum(1 for sp in sps if sp.lopsided)
    n_balanced = len(sps) - n_lopsided
    y_span = max(n_lopsided, 1) * pitch

    for sp in sps:
        fam_counts[sp.family] = fam_counts.get(sp.family, 0) + 1
        if sp.lopsided:
            angle = lop_i * angle_step
            x = math.cos(angle) * radius
            z = math.sin(angle) * radius
            y = lop_i * pitch
            r, g, b = FAM_COLORS[sp.family]
            helix += [x, y, z, r, g, b]
            lop_i += 1
        else:
            y = bal_i * (y_span / max(n_balanced, 1))
            r, g, b = FAM_COLORS["balanced"]
            axis += [0.0, y, 0.0, r, g, b]
            bal_i += 1

    n_total = lop_i + bal_i
    stats = {
        "total":         n_total,
        "lopsided":      lop_i,
        "balanced":      bal_i,
        "lopsided_pct":  100.0 * lop_i / max(n_total, 1),
        "balanced_pct":  100.0 * bal_i / max(n_total, 1),
        "families":      fam_counts,
    }
    return helix, axis, stats

# ── Shaders ───────────────────────────────────────────────────────────────────
VERT_SRC = """
#version 330 core
in vec3 position;
in vec3 color;
out vec3 v_color;
uniform mat4 mvp;
uniform float point_size;
void main() {
    gl_Position  = mvp * vec4(position, 1.0);
    gl_PointSize = point_size;
    v_color      = color;
}
"""

FRAG_SRC = """
#version 330 core
in  vec3 v_color;
out vec4 frag_color;
void main() {
    vec2 d = gl_PointCoord * 2.0 - 1.0;
    if (dot(d, d) > 1.0) discard;
    frag_color = vec4(v_color, 1.0);
}
"""

VERT_LINE = """
#version 330 core
in  vec3 position;
uniform mat4 mvp;
void main() { gl_Position = mvp * vec4(position, 1.0); }
"""
FRAG_LINE = """
#version 330 core
out vec4 frag_color;
uniform vec4 line_color;
void main() { frag_color = line_color; }
"""

# ── VBO helpers ───────────────────────────────────────────────────────────────

def _make_vbo(data: List[float]) -> Tuple[int, int]:
    """Upload float data to a VBO. Returns (vao_id, vbo_id)."""
    vao = gl.GLuint()
    gl.glGenVertexArrays(1, ctypes.byref(vao))
    gl.glBindVertexArray(vao)

    vbo = gl.GLuint()
    gl.glGenBuffers(1, ctypes.byref(vbo))
    gl.glBindBuffer(gl.GL_ARRAY_BUFFER, vbo)
    arr = (gl.GLfloat * len(data))(*data)
    gl.glBufferData(gl.GL_ARRAY_BUFFER, ctypes.sizeof(arr), arr, gl.GL_STATIC_DRAW)
    gl.glBindVertexArray(0)
    return vao.value, vbo.value


def _delete_vao_vbo(vao: int, vbo: int) -> None:
    gl.glDeleteVertexArrays(1, (gl.GLuint * 1)(vao))
    gl.glDeleteBuffers(1, (gl.GLuint * 1)(vbo))


def _bind_attrs(program_id: int, vao: int, stride: int = 6) -> None:
    """Bind position (loc 0) and color (loc 1) attributes, stride=6 floats."""
    gl.glBindVertexArray(vao)
    pos_loc   = gl.glGetAttribLocation(program_id, b"position")
    color_loc = gl.glGetAttribLocation(program_id, b"color")
    sz = ctypes.sizeof(gl.GLfloat)
    gl.glEnableVertexAttribArray(pos_loc)
    gl.glVertexAttribPointer(pos_loc,   3, gl.GL_FLOAT, False, stride * sz, 0)
    gl.glEnableVertexAttribArray(color_loc)
    gl.glVertexAttribPointer(color_loc, 3, gl.GL_FLOAT, False, stride * sz,
                             ctypes.c_void_p(3 * sz))
    gl.glBindVertexArray(0)

# ── Main window ───────────────────────────────────────────────────────────────

class HelixWindow(pyglet.window.Window):
    def __init__(self, all_sps: List[SP], args: argparse.Namespace):
        cfg = pyglet.gl.Config(
            double_buffer=True, depth_size=24,
            major_version=3, minor_version=3,
            forward_compatible=True,
        )
        super().__init__(
            width=1280, height=820,
            caption="PrimeHelix 3D Visualizer",
            resizable=True, config=cfg,
        )
        self.all_sps = all_sps
        self.visible_n  = min(args.n, len(all_sps))
        self.radius     = args.radius
        self.pitch      = args.pitch
        self.angle_step = args.angle_step
        self.rot_speed  = args.rotation_speed

        # View state
        self.rot_x     = -22.0
        self.rot_y     = 35.0
        self.zoom      = -230.0
        self.auto_rot  = True
        self._drag     = False
        self._last_t   = time.time()

        # Build shaders
        vert_shader  = Shader(VERT_SRC,  "vertex")
        frag_shader  = Shader(FRAG_SRC,  "fragment")
        vl_shader    = Shader(VERT_LINE, "vertex")
        fl_shader    = Shader(FRAG_LINE, "fragment")
        self._prog   = ShaderProgram(vert_shader, frag_shader)
        self._lprog  = ShaderProgram(vl_shader,   fl_shader)

        self._h_vao = self._h_vbo = 0
        self._a_vao = self._a_vbo = 0
        self._ax_vao = self._ax_vbo = 0
        self._stats: Dict = {}
        self._rebuild()

        # UI labels
        self._stat_lbl = pyglet.text.Label(
            "", x=12, y=self.height - 14, anchor_x="left", anchor_y="top",
            font_size=11, color=(220, 220, 255, 220),
            multiline=True, width=self.width - 24,
        )
        self._help_lbl = pyglet.text.Label(
            "drag=rotate  wheel=zoom  space=pause  [/] N  -/= radius  ,/. pitch  9/0 speed  R reset  S screenshot  Esc quit",
            x=12, y=14, anchor_x="left", anchor_y="bottom",
            font_size=10, color=(160, 160, 190, 180),
        )
        pyglet.clock.schedule_interval(self._tick, 1 / 60)

    # ── Data / geometry ───────────────────────────────────────────────────────

    def _rebuild(self) -> None:
        # Free old GPU resources
        for vao, vbo in ((self._h_vao, self._h_vbo),
                         (self._a_vao, self._a_vbo),
                         (self._ax_vao, self._ax_vbo)):
            if vao:
                _delete_vao_vbo(vao, vbo)

        sps = self.all_sps[: self.visible_n]
        helix, axis, stats = build_geometry(sps, self.radius, self.pitch, self.angle_step)
        self._stats = stats
        self._n_helix = len(helix) // 6
        self._n_axis  = len(axis)  // 6

        pid = self._prog.id

        if helix:
            self._h_vao, self._h_vbo = _make_vbo(helix)
            _bind_attrs(pid, self._h_vao)
        if axis:
            self._a_vao, self._a_vbo = _make_vbo(axis)
            _bind_attrs(pid, self._a_vao)

        # Y-axis guide line
        y_top = max(self._n_helix * self.pitch, 10.0)
        line_data = [0.0, 0.0, 0.0, 0.0, y_top * 1.08, 0.0]
        self._ax_vao, self._ax_vbo = _make_vbo(line_data)
        # Bind position only for line shader
        gl.glBindVertexArray(self._ax_vao)
        gl.glBindBuffer(gl.GL_ARRAY_BUFFER, self._ax_vbo)
        pos_loc = gl.glGetAttribLocation(self._lprog.id, b"position")
        sz = ctypes.sizeof(gl.GLfloat)
        gl.glEnableVertexAttribArray(pos_loc)
        gl.glVertexAttribPointer(pos_loc, 3, gl.GL_FLOAT, False, 3 * sz, 0)
        gl.glBindVertexArray(0)

        self._y_center = y_top / 2.0

    # ── MVP matrix ───────────────────────────────────────────────────────────

    def _mvp(self) -> pm.Mat4:
        w, h = self.get_size()
        proj = pm.Mat4.perspective_projection(
            aspect=w / max(h, 1), z_near=0.1, z_far=10_000.0, fov=55.0,
        )
        view = (
            pm.Mat4.from_translation(pm.Vec3(0.0, -self._y_center, self.zoom))
            @ pm.Mat4.from_rotation(math.radians(self.rot_x), pm.Vec3(1, 0, 0))
            @ pm.Mat4.from_rotation(math.radians(self.rot_y), pm.Vec3(0, 1, 0))
        )
        return proj @ view

    # ── Event handlers ───────────────────────────────────────────────────────

    def _tick(self, dt: float) -> None:
        if self.auto_rot:
            self.rot_y += self.rot_speed * dt * 35.0

    def on_draw(self) -> None:
        self.clear()
        gl.glClearColor(*BG)
        gl.glClear(gl.GL_COLOR_BUFFER_BIT | gl.GL_DEPTH_BUFFER_BIT)
        gl.glEnable(gl.GL_DEPTH_TEST)
        gl.glEnable(gl.GL_PROGRAM_POINT_SIZE)
        gl.glViewport(0, 0, *self.get_framebuffer_size())

        mvp = self._mvp()
        mvp_flat = (gl.GLfloat * 16)(*mvp)

        # Draw y-axis guide
        self._lprog.use()
        gl.glUniformMatrix4fv(
            gl.glGetUniformLocation(self._lprog.id, b"mvp"),
            1, False, mvp_flat,
        )
        gl.glUniform4f(
            gl.glGetUniformLocation(self._lprog.id, b"line_color"),
            0.35, 0.35, 0.40, 1.0,
        )
        if self._ax_vao:
            gl.glBindVertexArray(self._ax_vao)
            gl.glDrawArrays(gl.GL_LINES, 0, 2)
            gl.glBindVertexArray(0)
        self._lprog.stop()

        # Draw semiprime points
        self._prog.use()
        gl.glUniformMatrix4fv(
            gl.glGetUniformLocation(self._prog.id, b"mvp"),
            1, False, mvp_flat,
        )

        # Helix (lopsided)
        gl.glUniform1f(
            gl.glGetUniformLocation(self._prog.id, b"point_size"),
            DEFAULT_POINT_SIZE,
        )
        if self._h_vao and self._n_helix:
            gl.glBindVertexArray(self._h_vao)
            gl.glDrawArrays(gl.GL_POINTS, 0, self._n_helix)
            gl.glBindVertexArray(0)

        # Axis (balanced) — slightly larger
        gl.glUniform1f(
            gl.glGetUniformLocation(self._prog.id, b"point_size"),
            DEFAULT_POINT_SIZE * 1.8,
        )
        if self._a_vao and self._n_axis:
            gl.glBindVertexArray(self._a_vao)
            gl.glDrawArrays(gl.GL_POINTS, 0, self._n_axis)
            gl.glBindVertexArray(0)

        self._prog.stop()

        # 2-D overlay
        gl.glDisable(gl.GL_DEPTH_TEST)
        self._draw_ui()

    def _draw_ui(self) -> None:
        s = self._stats
        fam = s.get("families", {})
        self._stat_lbl.y     = self.height - 14
        self._stat_lbl.width = self.width - 24
        self._stat_lbl.text  = (
            f"N={s.get('total', 0):,} / {len(self.all_sps):,} computed  │  "
            f"lopsided {s.get('lopsided_pct', 0):.2f}%  balanced {s.get('balanced_pct', 0):.2f}%  │  "
            f"1x1 {fam.get('1x1', 0):,}  1x3 {fam.get('1x3', 0):,}  "
            f"3x3 {fam.get('3x3', 0):,}  even {fam.get('even', 0):,}  │  "
            f"radius {self.radius:.1f}  pitch {self.pitch:.3f}  speed {self.rot_speed:.2f}"
        )
        self._stat_lbl.draw()
        self._help_lbl.draw()

    def on_resize(self, width: int, height: int) -> bool:
        gl.glViewport(0, 0, *self.get_framebuffer_size())
        return pyglet.event.EVENT_HANDLED

    def on_mouse_drag(self, x, y, dx, dy, buttons, modifiers) -> None:
        if buttons & pyglet.window.mouse.LEFT:
            self.rot_y += dx * 0.35
            self.rot_x -= dy * 0.35

    def on_mouse_scroll(self, x, y, sx, sy) -> None:
        self.zoom = min(-20.0, max(-4_000.0, self.zoom + sy * 10.0))

    def on_key_press(self, symbol, modifiers) -> None:
        K = pyglet.window.key
        step = max(200, self.visible_n // 10)
        if symbol == K.ESCAPE:
            self.close()
        elif symbol == K.SPACE:
            self.auto_rot = not self.auto_rot
        elif symbol == K.BRACKETLEFT:
            self.visible_n = max(10, self.visible_n - step)
            self._rebuild()
        elif symbol == K.BRACKETRIGHT:
            self.visible_n = min(len(self.all_sps), self.visible_n + step)
            self._rebuild()
        elif symbol == K.MINUS:
            self.radius = max(5.0, self.radius - 5.0)
            self._rebuild()
        elif symbol == K.EQUAL:
            self.radius += 5.0
            self._rebuild()
        elif symbol == K.COMMA:
            self.pitch = max(0.01, self.pitch - 0.02)
            self._rebuild()
        elif symbol == K.PERIOD:
            self.pitch += 0.02
            self._rebuild()
        elif symbol == K._9:
            self.rot_speed = max(0.0, self.rot_speed - 0.05)
        elif symbol == K._0:
            self.rot_speed += 0.05
        elif symbol == K.R:
            self.rot_x, self.rot_y, self.zoom = -22.0, 35.0, -230.0
        elif symbol == K.S:
            fname = f"primehelix_{int(time.time())}.png"
            pyglet.image.get_buffer_manager().get_color_buffer().save(fname)
            print(f"screenshot saved: {fname}")

# ── CLI ───────────────────────────────────────────────────────────────────────

def _parse() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="PrimeHelix 3D native visualizer")
    p.add_argument("--limit",          type=int,   default=DEFAULT_LIMIT,
                   help="Compute semiprimes n=p·q up to this value (default 1M)")
    p.add_argument("--n",              type=int,   default=DEFAULT_N,
                   help="Initial number of semiprimes to display")
    p.add_argument("--bit-gap",        type=int,   default=8,
                   help="Lopsided threshold: floor(log2 q)−floor(log2 p) > bit-gap")
    p.add_argument("--radius",         type=float, default=DEFAULT_RADIUS)
    p.add_argument("--pitch",          type=float, default=DEFAULT_PITCH)
    p.add_argument("--angle-step",     type=float, default=DEFAULT_ANGLE_STEP)
    p.add_argument("--rotation-speed", type=float, default=DEFAULT_SPEED)
    return p.parse_args()


def main() -> None:
    args = _parse()
    print(f"Generating semiprimes ≤ {args.limit:,} …")
    sps = generate_semiprimes(args.limit, bit_gap=args.bit_gap)
    if not sps:
        sys.exit("No semiprimes generated — try a larger --limit.")
    print(f"Generated {len(sps):,} semiprimes.")
    win = HelixWindow(sps, args)
    pyglet.app.run()


if __name__ == "__main__":
    main()
