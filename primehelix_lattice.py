#!/usr/bin/env python3
"""
primehelix_lattice.py — 3D semiprime lattice visualizer with family edges

Install:  pip install pyglet
Run:      python primehelix_lattice.py
          python primehelix_lattice.py --n 600 --radius 80 --speed 0.3

Controls:
  Left drag      rotate          Scroll       zoom
  Double-click   reset view      Space        toggle auto-rotate
  1 / 2 / 3      toggle family strands / axis spokes / cross chords
  [ / ]          N −100 / +100
  R              reset camera    S            screenshot    Esc  quit
"""
from __future__ import annotations

import argparse
import ctypes
import math
import sys
import time
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

try:
    import pyglet
    from pyglet import gl
    from pyglet.graphics.shader import Shader, ShaderProgram
    import pyglet.math as pm
except ImportError:
    sys.exit("pyglet not found — run: pip install pyglet")

# ── Palette ───────────────────────────────────────────────────────────────────
BG = (10 / 255, 10 / 255, 18 / 255, 1.0)

FAM_RGB: Dict[str, Tuple[float, float, float]] = {
    "1x1":      (0.55, 0.20, 0.85),
    "1x3":      (0.20, 0.82, 0.35),
    "3x3":      (1.00, 0.42, 0.28),
    "even":     (0.62, 0.62, 0.62),
    "balanced": (1.00, 0.75, 0.10),
}

STRAND_ALPHA = 0.70
SPOKE_ALPHA  = 0.50
CHORD_ALPHA  = 0.30
CHORD_RGB    = (0.70, 0.70, 0.75)

POINT_SIZE_HELIX = 4.0
POINT_SIZE_AXIS  = 8.0

# ── Number theory (pure Python) ───────────────────────────────────────────────

def _is_prime(n: int) -> bool:
    if n < 2: return False
    if n < 4: return True
    if n % 2 == 0 or n % 3 == 0: return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0: return False
        i += 6
    return True


def _spf(n: int) -> int:
    """Smallest prime factor; returns n if n is prime."""
    if n % 2 == 0: return 2
    i = 3
    while i * i <= n:
        if n % i == 0: return i
        i += 2
    return n


def _residue_family(p: int, q: int) -> str:
    if p == 2 or q == 2: return "even"
    a, b = p % 4, q % 4
    if a == 1 and b == 1: return "1x1"
    if a == 3 and b == 3: return "3x3"
    return "1x3"


def _lopsided(p: int, q: int) -> bool:
    return (q.bit_length() - p.bit_length()) > 8

# ── Semiprime dataclass ───────────────────────────────────────────────────────

@dataclass
class SP:
    idx:     int            # position in full N-sorted list
    n:       int
    p:       int
    q:       int
    lopsided: bool
    family:  str
    x: float = 0.0
    y: float = 0.0
    z: float = 0.0
    lop_idx: int = -1       # index among lopsided semiprimes
    bal_idx: int = -1       # index among balanced semiprimes

# ── Data generation ───────────────────────────────────────────────────────────

def generate_semiprimes(count: int) -> List[SP]:
    """Scan integers from 4 upward; return first `count` semiprimes."""
    out: List[SP] = []
    k = 4
    while len(out) < count:
        s = _spf(k)
        if s < k:
            q = k // s
            if _is_prime(q):
                out.append(SP(
                    idx=len(out), n=k, p=s, q=q,
                    lopsided=_lopsided(s, q),
                    family=_residue_family(s, q),
                ))
        k += 1
    return out

# ── Geometry ──────────────────────────────────────────────────────────────────

def assign_positions(sps: List[SP], radius: float, span: float) -> None:
    """Assign (x, y, z) to each SP in-place and set lop_idx / bal_idx."""
    lop = [s for s in sps if s.lopsided]
    bal = [s for s in sps if not s.lopsided]

    for i, sp in enumerate(lop):
        sp.lop_idx = i
        angle = i * 0.27
        sp.x = math.cos(angle) * radius
        sp.z = math.sin(angle) * radius
        sp.y = (i / max(len(lop) - 1, 1) - 0.5) * span

    for i, sp in enumerate(bal):
        sp.bal_idx = i
        sp.x = 0.0
        sp.z = 0.0
        sp.y = (i / max(len(bal) - 1, 1) - 0.5) * span

# ── Edge builders ─────────────────────────────────────────────────────────────

def _vert6(sp: SP, r: float, g: float, b: float, a: float) -> List[float]:
    return [sp.x, sp.y, sp.z, r, g, b, a]


def build_strands(sps: List[SP]) -> List[float]:
    """Family strands: line between consecutive semiprimes of the same family."""
    verts: List[float] = []
    prev: Dict[str, SP] = {}
    for sp in sps:
        r, g, b = FAM_RGB[sp.family if sp.lopsided else "balanced"]
        if sp.family in prev:
            p = prev[sp.family]
            pr, pg, pb = FAM_RGB[p.family if p.lopsided else "balanced"]
            verts += _vert6(p, pr, pg, pb, STRAND_ALPHA)
            verts += _vert6(sp, r, g, b, STRAND_ALPHA)
        prev[sp.family] = sp
    return verts


def build_spokes(sps: List[SP]) -> List[float]:
    """Axis spokes: balanced → nearest lopsided by index distance."""
    lop = [s for s in sps if s.lopsided]
    if not lop:
        return []
    r, g, b = FAM_RGB["balanced"]
    verts: List[float] = []
    for sp in sps:
        if sp.lopsided:
            continue
        # find nearest lopsided by list index
        nearest: Optional[SP] = None
        best = 10 ** 9
        for lp in lop:
            d = abs(lp.idx - sp.idx)
            if d < best:
                best = d
                nearest = lp
        if nearest is None:
            continue
        lr, lg, lb = FAM_RGB[nearest.family]
        verts += _vert6(sp,      r,  g,  b,  SPOKE_ALPHA)
        verts += _vert6(nearest, lr, lg, lb, SPOKE_ALPHA)
    return verts


def build_chords(sps: List[SP]) -> List[float]:
    """Cross chords: adjacent lopsided semiprimes of different families."""
    lop = [s for s in sps if s.lopsided]
    verts: List[float] = []
    cr, cg, cb = CHORD_RGB
    for a, b_sp in zip(lop, lop[1:]):
        if a.family != b_sp.family:
            verts += _vert6(a,    cr, cg, cb, CHORD_ALPHA)
            verts += _vert6(b_sp, cr, cg, cb, CHORD_ALPHA)
    return verts


def compute_stats(sps: List[SP]) -> Dict:
    n = len(sps)
    lop = sum(1 for s in sps if s.lopsided)
    fam: Dict[str, int] = {}
    for s in sps:
        fam[s.family] = fam.get(s.family, 0) + 1
    return {
        "n": n, "lopsided": lop, "balanced": n - lop,
        "lop_pct": 100.0 * lop / max(n, 1),
        "bal_pct": 100.0 * (n - lop) / max(n, 1),
        "fam": fam,
    }

# ── Shaders ───────────────────────────────────────────────────────────────────

_PT_VERT = """
#version 330 core
in vec3 position;
in vec3 color;
out vec3 v_color;
uniform mat4 mvp;
uniform float pt_size;
void main() {
    gl_Position  = mvp * vec4(position, 1.0);
    gl_PointSize = pt_size;
    v_color      = color;
}
"""

_PT_FRAG = """
#version 330 core
in  vec3 v_color;
out vec4 frag_color;
void main() {
    vec2 d = gl_PointCoord * 2.0 - 1.0;
    if (dot(d, d) > 1.0) discard;
    frag_color = vec4(v_color, 1.0);
}
"""

_LN_VERT = """
#version 330 core
in vec3 position;
in vec4 color;
out vec4 v_color;
uniform mat4 mvp;
void main() {
    gl_Position = mvp * vec4(position, 1.0);
    v_color     = color;
}
"""

_LN_FRAG = """
#version 330 core
in  vec4 v_color;
out vec4 frag_color;
void main() { frag_color = v_color; }
"""

# ── VBO utilities ─────────────────────────────────────────────────────────────

class VBO:
    __slots__ = ("vao", "vbo", "count", "stride")

    def __init__(self, data: List[float], stride: int):
        self.stride = stride
        self.count = len(data) // stride
        vao = gl.GLuint()
        vbo = gl.GLuint()
        gl.glGenVertexArrays(1, ctypes.byref(vao))
        gl.glGenBuffers(1, ctypes.byref(vbo))
        self.vao = vao.value
        self.vbo = vbo.value
        if data:
            arr = (gl.GLfloat * len(data))(*data)
            gl.glBindBuffer(gl.GL_ARRAY_BUFFER, self.vbo)
            gl.glBufferData(gl.GL_ARRAY_BUFFER, ctypes.sizeof(arr), arr, gl.GL_STATIC_DRAW)
            gl.glBindBuffer(gl.GL_ARRAY_BUFFER, 0)

    def delete(self) -> None:
        gl.glDeleteVertexArrays(1, (gl.GLuint * 1)(self.vao))
        gl.glDeleteBuffers(1, (gl.GLuint * 1)(self.vbo))

    def bind_point_attrs(self, prog_id: int) -> None:
        """Bind position(vec3) + color(vec3), stride 6 floats."""
        sz = ctypes.sizeof(gl.GLfloat)
        gl.glBindVertexArray(self.vao)
        gl.glBindBuffer(gl.GL_ARRAY_BUFFER, self.vbo)
        for name, n_comp, offset in (
            (b"position", 3, 0),
            (b"color",    3, 3 * sz),
        ):
            loc = gl.glGetAttribLocation(prog_id, name)
            if loc >= 0:
                gl.glEnableVertexAttribArray(loc)
                gl.glVertexAttribPointer(loc, n_comp, gl.GL_FLOAT, False,
                                         self.stride * sz,
                                         ctypes.c_void_p(offset))
        gl.glBindVertexArray(0)

    def bind_line_attrs(self, prog_id: int) -> None:
        """Bind position(vec3) + color(vec4), stride 7 floats."""
        sz = ctypes.sizeof(gl.GLfloat)
        gl.glBindVertexArray(self.vao)
        gl.glBindBuffer(gl.GL_ARRAY_BUFFER, self.vbo)
        for name, n_comp, offset in (
            (b"position", 3, 0),
            (b"color",    4, 3 * sz),
        ):
            loc = gl.glGetAttribLocation(prog_id, name)
            if loc >= 0:
                gl.glEnableVertexAttribArray(loc)
                gl.glVertexAttribPointer(loc, n_comp, gl.GL_FLOAT, False,
                                         self.stride * sz,
                                         ctypes.c_void_p(offset))
        gl.glBindVertexArray(0)


def _upload_points(sps: List[SP], filt_lopsided: Optional[bool]) -> List[float]:
    data: List[float] = []
    for sp in sps:
        if filt_lopsided is not None and sp.lopsided != filt_lopsided:
            continue
        if sp.lopsided:
            r, g, b = FAM_RGB[sp.family]
        else:
            r, g, b = FAM_RGB["balanced"]
        data += [sp.x, sp.y, sp.z, r, g, b]
    return data

# ── MVP helper ────────────────────────────────────────────────────────────────

def _mvp(rot_x: float, rot_y: float, zoom: float,
         y_center: float, w: int, h: int) -> pm.Mat4:
    proj = pm.Mat4.perspective_projection(
        aspect=w / max(h, 1), z_near=0.1, z_far=50_000.0, fov=55.0,
    )
    view = (
        pm.Mat4.from_translation(pm.Vec3(0.0, -y_center, zoom))
        @ pm.Mat4.from_rotation(math.radians(rot_x), pm.Vec3(1, 0, 0))
        @ pm.Mat4.from_rotation(math.radians(rot_y), pm.Vec3(0, 1, 0))
    )
    return proj @ view

# ── Window ────────────────────────────────────────────────────────────────────

class LatticeWindow(pyglet.window.Window):

    _RESET_ROT_X  = -20.0
    _RESET_ROT_Y  =  30.0
    _RESET_ZOOM   = -300.0
    _DOUBLE_CLICK = 0.35   # seconds

    def __init__(self, args: argparse.Namespace):
        cfg = pyglet.gl.Config(
            double_buffer=True, depth_size=24,
            major_version=3, minor_version=3,
            forward_compatible=True,
        )
        super().__init__(
            width=1280, height=820,
            caption="PrimeHelix Lattice",
            resizable=True, config=cfg,
        )
        self._args    = args
        self._count   = args.n
        self._radius  = args.radius
        self._speed   = args.speed

        self._rot_x   = self._RESET_ROT_X
        self._rot_y   = self._RESET_ROT_Y
        self._zoom    = self._RESET_ZOOM
        self._auto    = True
        self._last_t  = time.time()
        self._click_t = 0.0

        # Edge toggles
        self._show_strands = True
        self._show_spokes  = True
        self._show_chords  = True

        # GPU handles (populated in _rebuild)
        self._pt_prog: Optional[ShaderProgram] = None
        self._ln_prog: Optional[ShaderProgram] = None
        self._vbo_helix:   Optional[VBO] = None
        self._vbo_axis:    Optional[VBO] = None
        self._vbo_strands: Optional[VBO] = None
        self._vbo_spokes:  Optional[VBO] = None
        self._vbo_chords:  Optional[VBO] = None
        self._stats: Dict = {}
        self._y_center = 0.0

        self._build_shaders()
        self._rebuild()
        self._build_ui()
        pyglet.clock.schedule_interval(self._tick, 1 / 60)

    # ── Build ─────────────────────────────────────────────────────────────────

    def _build_shaders(self) -> None:
        self._pt_prog = ShaderProgram(Shader(_PT_VERT, "vertex"),
                                      Shader(_PT_FRAG, "fragment"))
        self._ln_prog = ShaderProgram(Shader(_LN_VERT, "vertex"),
                                      Shader(_LN_FRAG, "fragment"))

    def _free_vbos(self) -> None:
        for attr in ("_vbo_helix", "_vbo_axis",
                     "_vbo_strands", "_vbo_spokes", "_vbo_chords"):
            v = getattr(self, attr)
            if v is not None:
                v.delete()
                setattr(self, attr, None)

    def _rebuild(self) -> None:
        self._free_vbos()
        sps = generate_semiprimes(self._count)
        span = self._radius * 2.8
        assign_positions(sps, self._radius, span)
        self._stats = compute_stats(sps)

        y_vals = [s.y for s in sps]
        self._y_center = (max(y_vals) + min(y_vals)) / 2 if y_vals else 0.0

        pid_pt = self._pt_prog.id
        pid_ln = self._ln_prog.id

        # Points
        h_data = _upload_points(sps, filt_lopsided=True)
        a_data = _upload_points(sps, filt_lopsided=False)
        self._vbo_helix = VBO(h_data, stride=6)
        self._vbo_axis  = VBO(a_data, stride=6)
        if h_data: self._vbo_helix.bind_point_attrs(pid_pt)
        if a_data: self._vbo_axis.bind_point_attrs(pid_pt)

        # Edges
        s_data = build_strands(sps)
        sp_data = build_spokes(sps)
        c_data  = build_chords(sps)
        self._vbo_strands = VBO(s_data,  stride=7)
        self._vbo_spokes  = VBO(sp_data, stride=7)
        self._vbo_chords  = VBO(c_data,  stride=7)
        for vbo in (self._vbo_strands, self._vbo_spokes, self._vbo_chords):
            if vbo.count:
                vbo.bind_line_attrs(pid_ln)

    def _build_ui(self) -> None:
        self._stat_lbl = pyglet.text.Label(
            "", x=12, y=14, anchor_x="left", anchor_y="bottom",
            font_size=11, color=(210, 210, 240, 220),
            multiline=True, width=self.width - 24,
        )
        self._help_lbl = pyglet.text.Label(
            "drag=rotate  scroll=zoom  dbl-click=reset  space=auto  "
            "1/2/3=edges  [/]=N  R=camera  S=screenshot",
            x=self.width - 12, y=14, anchor_x="right", anchor_y="bottom",
            font_size=10, color=(130, 130, 160, 170),
        )

    # ── Tick / draw ───────────────────────────────────────────────────────────

    def _tick(self, dt: float) -> None:
        if self._auto:
            self._rot_y += self._speed * dt * 35.0

    def on_draw(self) -> None:
        self.clear()
        gl.glClearColor(*BG)
        gl.glClear(gl.GL_COLOR_BUFFER_BIT | gl.GL_DEPTH_BUFFER_BIT)
        gl.glEnable(gl.GL_DEPTH_TEST)
        gl.glEnable(gl.GL_PROGRAM_POINT_SIZE)
        gl.glEnable(gl.GL_BLEND)
        gl.glBlendFunc(gl.GL_SRC_ALPHA, gl.GL_ONE_MINUS_SRC_ALPHA)
        gl.glViewport(0, 0, *self.get_framebuffer_size())

        w, h = self.get_size()
        mvp = _mvp(self._rot_x, self._rot_y, self._zoom, self._y_center, w, h)
        mvp_flat = (gl.GLfloat * 16)(*mvp)

        # ── Lines ─────────────────────────────────────────────────────────────
        self._ln_prog.use()
        mvp_loc = gl.glGetUniformLocation(self._ln_prog.id, b"mvp")
        gl.glUniformMatrix4fv(mvp_loc, 1, False, mvp_flat)

        pairs = [
            (self._vbo_strands, self._show_strands),
            (self._vbo_spokes,  self._show_spokes),
            (self._vbo_chords,  self._show_chords),
        ]
        for vbo, visible in pairs:
            if visible and vbo and vbo.count:
                gl.glBindVertexArray(vbo.vao)
                gl.glDrawArrays(gl.GL_LINES, 0, vbo.count)
                gl.glBindVertexArray(0)
        self._ln_prog.stop()

        # ── Points ────────────────────────────────────────────────────────────
        self._pt_prog.use()
        mvp_loc_pt = gl.glGetUniformLocation(self._pt_prog.id, b"mvp")
        sz_loc     = gl.glGetUniformLocation(self._pt_prog.id, b"pt_size")
        gl.glUniformMatrix4fv(mvp_loc_pt, 1, False, mvp_flat)

        gl.glUniform1f(sz_loc, POINT_SIZE_HELIX)
        if self._vbo_helix and self._vbo_helix.count:
            gl.glBindVertexArray(self._vbo_helix.vao)
            gl.glDrawArrays(gl.GL_POINTS, 0, self._vbo_helix.count)
            gl.glBindVertexArray(0)

        gl.glUniform1f(sz_loc, POINT_SIZE_AXIS)
        if self._vbo_axis and self._vbo_axis.count:
            gl.glBindVertexArray(self._vbo_axis.vao)
            gl.glDrawArrays(gl.GL_POINTS, 0, self._vbo_axis.count)
            gl.glBindVertexArray(0)

        self._pt_prog.stop()

        # ── 2-D overlay ───────────────────────────────────────────────────────
        gl.glDisable(gl.GL_DEPTH_TEST)
        self._update_labels()
        self._stat_lbl.draw()
        self._help_lbl.draw()

    def _update_labels(self) -> None:
        s = self._stats
        fam = s.get("fam", {})
        e = (
            f"[1]strands={'on' if self._show_strands else 'off'}  "
            f"[2]spokes={'on' if self._show_spokes else 'off'}  "
            f"[3]chords={'on' if self._show_chords else 'off'}"
        )
        self._stat_lbl.text = (
            f"N={s.get('n', 0):,}  lopsided {s.get('lop_pct', 0):.1f}%  "
            f"balanced {s.get('bal_pct', 0):.1f}%  │  "
            f"1x1={fam.get('1x1', 0):,}  1x3={fam.get('1x3', 0):,}  "
            f"3x3={fam.get('3x3', 0):,}  even={fam.get('even', 0):,}  │  "
            f"radius={self._radius:.0f}  zoom={abs(self._zoom):.0f}  {e}"
        )
        self._help_lbl.x = self.width - 12

    # ── Events ────────────────────────────────────────────────────────────────

    def on_resize(self, w: int, h: int) -> bool:
        gl.glViewport(0, 0, *self.get_framebuffer_size())
        self._help_lbl.x = w - 12
        self._stat_lbl.width = w - 24
        return pyglet.event.EVENT_HANDLED

    def on_mouse_press(self, x, y, button, modifiers) -> None:
        if button != pyglet.window.mouse.LEFT:
            return
        now = time.time()
        if now - self._click_t < self._DOUBLE_CLICK:
            self._reset_view()
        self._click_t = now

    def on_mouse_drag(self, x, y, dx, dy, buttons, modifiers) -> None:
        if buttons & pyglet.window.mouse.LEFT:
            self._rot_y += dx * 0.35
            self._rot_x -= dy * 0.35

    def on_mouse_scroll(self, x, y, sx, sy) -> None:
        self._zoom = min(-20.0, max(-8_000.0, self._zoom + sy * 15.0))

    def on_key_press(self, symbol, modifiers) -> None:
        K = pyglet.window.key
        if symbol == K.ESCAPE:
            self.close()
        elif symbol == K.SPACE:
            self._auto = not self._auto
        elif symbol == K._1:
            self._show_strands = not self._show_strands
        elif symbol == K._2:
            self._show_spokes = not self._show_spokes
        elif symbol == K._3:
            self._show_chords = not self._show_chords
        elif symbol == K.BRACKETLEFT:
            self._count = max(10, self._count - 100)
            self._rebuild()
        elif symbol == K.BRACKETRIGHT:
            self._count += 100
            self._rebuild()
        elif symbol == K.R:
            self._reset_view()
        elif symbol == K.S:
            fname = f"lattice_{int(time.time())}.png"
            pyglet.image.get_buffer_manager().get_color_buffer().save(fname)
            print(f"screenshot: {fname}")

    def _reset_view(self) -> None:
        self._rot_x = self._RESET_ROT_X
        self._rot_y = self._RESET_ROT_Y
        self._zoom  = self._RESET_ZOOM

# ── CLI / main ────────────────────────────────────────────────────────────────

def _parse() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="PrimeHelix 3D lattice visualizer")
    p.add_argument("--n",      type=int,   default=600,
                   help="Initial semiprime count (use ≥2000 for a visible helix)")
    p.add_argument("--radius", type=float, default=80.0, help="Helix radius")
    p.add_argument("--speed",  type=float, default=0.3,  help="Auto-rotation speed")
    return p.parse_args()


def main() -> None:
    args = _parse()
    print(f"Building {args.n} semiprimes …")
    t0 = time.time()
    sps = generate_semiprimes(args.n)
    lop = sum(1 for s in sps if s.lopsided)
    lop_pct = 100.0 * lop / max(len(sps), 1)
    print(f"Done in {time.time()-t0:.2f}s — "
          f"{len(sps)} semiprimes, {lop_pct:.1f}% lopsided")
    if lop < 20:
        print(f"  ⚠  Only {lop} lopsided semiprimes visible at N={args.n}.")
        print(f"     Run with --n 2000 or higher for a populated helix.")
    win = LatticeWindow(args)
    pyglet.app.run()


if __name__ == "__main__":
    main()
