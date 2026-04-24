#!/usr/bin/env python3
"""
primehelix_lattice.py — 3D semiprime lattice visualizer with family edges

Install:  pip install pyglet
Run:      python primehelix_lattice.py
          python primehelix_lattice.py --n 2000 --radius 80 --speed 0.3

Controls:
  Left drag      rotate          Scroll       zoom in/out
  Double-click   reset view      Space        toggle auto-rotate
  1 / 2 / 3      toggle family strands / axis spokes / cross chords
  [ / ]          N −100 / +100  R            reset camera
  S              screenshot      Esc          quit
"""
from __future__ import annotations

import array
import argparse
import bisect
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

# ── Constants ─────────────────────────────────────────────────────────────────
# Sieve 1M integers → ~210k semiprimes at ~73% lopsided (bit-gap criterion).
SCAN_LIMIT = 1_000_000

BG = (10 / 255, 10 / 255, 18 / 255, 1.0)  # #0a0a12

FAM_RGB: Dict[str, Tuple[float, float, float]] = {
    "1x1":      (159 / 255, 143 / 255, 248 / 255),  # #9f8ff8  purple
    "1x3":      ( 62 / 255, 207 / 255, 150 / 255),  # #3ecf96  green
    "3x3":      (240 / 255, 112 / 255,  80 / 255),  # #f07050  coral
    "even":     (170 / 255, 170 / 255, 170 / 255),  # #aaaaaa  gray
    "balanced": (255 / 255, 208 / 255, 128 / 255),  # #ffd080  amber
}

STRAND_ALPHA = 0.70
SPOKE_ALPHA  = 0.50
CHORD_ALPHA  = 0.30
CHORD_RGB    = (0.70, 0.70, 0.75)

POINT_SIZE_HELIX = 5.5
POINT_SIZE_AXIS  = 9.5

# ── Number theory ─────────────────────────────────────────────────────────────

def _lopsided(p: int, q: int) -> bool:
    # balance = (q-p)/sqrt(pq) >= 10  ↔  (q-p)² >= 100·p·q  (integer, no sqrt)
    d = q - p
    return d * d >= 100 * p * q


def _residue_family(p: int, q: int) -> str:
    if p == 2 or q == 2:
        return "even"
    a, b = p % 4, q % 4
    if a == 1 and b == 1:
        return "1x1"
    if a == 3 and b == 3:
        return "3x3"
    return "1x3"

# ── Semiprime dataclass ───────────────────────────────────────────────────────

@dataclass
class SP:
    idx:      int
    n:        int
    p:        int
    q:        int
    lopsided: bool
    family:   str
    x: float = 0.0
    y: float = 0.0
    z: float = 0.0
    lop_idx: int = -1
    bal_idx: int = -1

# ── Sieve-based semiprime generation ─────────────────────────────────────────
# The fix for the low-lopsided-fraction bug: instead of scanning from n=4
# and stopping at count semiprimes (which puts us in the n<11k range where
# bit-gap lopsided fraction is ~25%), we sieve the full 1M range once and
# uniformly sample `count` from it. The 1M range has ~73% lopsided fraction.

_ALL_SP_DATA: Optional[List[Tuple[int, int, int, bool, str]]] = None


def _build_sp_data(limit: int) -> List[Tuple[int, int, int, bool, str]]:
    """Smallest-prime-factor sieve; return all semiprime tuples in [4, limit]."""
    spf = array.array('i', range(limit + 1))
    i = 2
    while i * i <= limit:
        if spf[i] == i:
            for j in range(i * i, limit + 1, i):
                if spf[j] == j:
                    spf[j] = i
        i += 1
    out: List[Tuple[int, int, int, bool, str]] = []
    for n in range(4, limit + 1):
        p = spf[n]
        if p == n:
            continue          # n is prime
        q = n // p
        if spf[q] != q:
            continue          # q is not prime → not a semiprime
        out.append((n, p, q, _lopsided(p, q), _residue_family(p, q)))
    return out


def generate_semiprimes(count: int) -> List[SP]:
    """Return `count` semiprimes sampled uniformly from [4, SCAN_LIMIT].

    Uniform sampling from 1M integers preserves the ~73% lopsided fraction
    regardless of how small `count` is. The sieve is built once and cached.
    """
    global _ALL_SP_DATA
    if _ALL_SP_DATA is None:
        print(f"Sieving to {SCAN_LIMIT:,} …", end=" ", flush=True)
        t0 = time.time()
        _ALL_SP_DATA = _build_sp_data(SCAN_LIMIT)
        print(f"{len(_ALL_SP_DATA):,} semiprimes ({time.time() - t0:.2f}s)")

    data = _ALL_SP_DATA
    if len(data) > count:
        step = len(data) / count
        data = [data[int(i * step)] for i in range(count)]

    return [SP(i, n, p, q, lop, fam)
            for i, (n, p, q, lop, fam) in enumerate(data)]

# ── Geometry ──────────────────────────────────────────────────────────────────

def assign_positions(sps: List[SP], radius: float, span: float) -> None:
    lop = [s for s in sps if s.lopsided]
    bal = [s for s in sps if not s.lopsided]
    n_lop = len(lop)
    n_bal = len(bal)

    for i, sp in enumerate(lop):
        sp.lop_idx = i
        angle = i * 0.27
        sp.x = math.cos(angle) * radius
        sp.z = math.sin(angle) * radius
        sp.y = (i / max(n_lop - 1, 1) - 0.5) * span

    for i, sp in enumerate(bal):
        sp.bal_idx = i
        sp.x = 0.0
        sp.z = 0.0
        sp.y = (i / max(n_bal - 1, 1) - 0.5) * span

# ── Edge / line builders ──────────────────────────────────────────────────────

def _vert7(sp: SP, r: float, g: float, b: float, a: float) -> List[float]:
    return [sp.x, sp.y, sp.z, r, g, b, a]


def build_strands(sps: List[SP]) -> List[float]:
    """Colored lines between consecutive lopsided semiprimes of the same family."""
    verts: List[float] = []
    prev: Dict[str, SP] = {}
    for sp in sps:
        if not sp.lopsided:
            continue
        r, g, b = FAM_RGB[sp.family]
        if sp.family in prev:
            p = prev[sp.family]
            pr, pg, pb = FAM_RGB[p.family]
            verts += _vert7(p,  pr, pg, pb, STRAND_ALPHA)
            verts += _vert7(sp, r,  g,  b,  STRAND_ALPHA)
        prev[sp.family] = sp
    return verts


def build_spokes(sps: List[SP]) -> List[float]:
    """Amber lines from each balanced point to the nearest lopsided by index."""
    lop = [s for s in sps if s.lopsided]
    if not lop:
        return []
    lop_indices = [s.idx for s in lop]   # already sorted ascending

    r, g, b = FAM_RGB["balanced"]
    verts: List[float] = []
    for sp in sps:
        if sp.lopsided:
            continue
        pos = bisect.bisect_left(lop_indices, sp.idx)
        candidates: List[SP] = []
        if pos < len(lop):
            candidates.append(lop[pos])
        if pos > 0:
            candidates.append(lop[pos - 1])
        if not candidates:
            continue
        nearest = min(candidates, key=lambda lp: abs(lp.idx - sp.idx))
        lr, lg, lb = FAM_RGB[nearest.family]
        verts += _vert7(sp,      r,  g,  b,  SPOKE_ALPHA)
        verts += _vert7(nearest, lr, lg, lb, SPOKE_ALPHA)
    return verts


def build_chords(sps: List[SP]) -> List[float]:
    """Gray lines between adjacent lopsided semiprimes of different families."""
    lop = [s for s in sps if s.lopsided]
    verts: List[float] = []
    cr, cg, cb = CHORD_RGB
    for a_sp, b_sp in zip(lop, lop[1:]):
        if a_sp.family != b_sp.family:
            verts += _vert7(a_sp, cr, cg, cb, CHORD_ALPHA)
            verts += _vert7(b_sp, cr, cg, cb, CHORD_ALPHA)
    return verts


def build_axis_line(span: float) -> List[float]:
    """Thin amber vertical line along the balanced-point axis."""
    r, g, b = FAM_RGB["balanced"]
    a = 0.45
    return [
        0.0, -span / 2, 0.0, r, g, b, a,
        0.0,  span / 2, 0.0, r, g, b, a,
    ]


def compute_stats(sps: List[SP]) -> Dict:
    total = len(sps)
    lop = sum(1 for s in sps if s.lopsided)
    fam: Dict[str, int] = {}
    for s in sps:
        fam[s.family] = fam.get(s.family, 0) + 1
    return {
        "n": total,
        "lopsided": lop,
        "balanced": total - lop,
        "lop_pct": 100.0 * lop / max(total, 1),
        "bal_pct": 100.0 * (total - lop) / max(total, 1),
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

# ── VBO helpers ───────────────────────────────────────────────────────────────

class VBO:
    __slots__ = ("vao", "vbo", "count", "stride")

    def __init__(self, data: List[float], stride: int):
        self.stride = stride
        self.count  = len(data) // stride
        vao = gl.GLuint()
        vbo = gl.GLuint()
        gl.glGenVertexArrays(1, ctypes.byref(vao))
        gl.glGenBuffers(1, ctypes.byref(vbo))
        self.vao = vao.value
        self.vbo = vbo.value
        if data:
            arr = (gl.GLfloat * len(data))(*data)
            gl.glBindBuffer(gl.GL_ARRAY_BUFFER, self.vbo)
            gl.glBufferData(gl.GL_ARRAY_BUFFER, ctypes.sizeof(arr),
                            arr, gl.GL_STATIC_DRAW)
            gl.glBindBuffer(gl.GL_ARRAY_BUFFER, 0)

    def delete(self) -> None:
        gl.glDeleteVertexArrays(1, (gl.GLuint * 1)(self.vao))
        gl.glDeleteBuffers(1, (gl.GLuint * 1)(self.vbo))

    def bind_point_attrs(self, prog_id: int) -> None:
        """Bind position(vec3) + color(vec3), stride=6."""
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
        """Bind position(vec3) + color(vec4), stride=7."""
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
    _DOUBLE_CLICK = 0.35

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
        self._args   = args
        self._count  = args.n
        self._radius = args.radius
        self._speed  = args.speed

        # Default zoom fits the full helix in view: span = radius*4, FOV≈55°
        self._zoom_default = -(args.radius * 4.0)
        self._rot_x  = self._RESET_ROT_X
        self._rot_y  = self._RESET_ROT_Y
        self._zoom   = self._zoom_default
        self._auto   = True
        self._click_t = 0.0

        self._show_strands = True
        self._show_spokes  = True
        self._show_chords  = True

        self._pt_prog:      Optional[ShaderProgram] = None
        self._ln_prog:      Optional[ShaderProgram] = None
        self._vbo_helix:    Optional[VBO] = None
        self._vbo_axis_pts: Optional[VBO] = None
        self._vbo_axis_ln:  Optional[VBO] = None
        self._vbo_strands:  Optional[VBO] = None
        self._vbo_spokes:   Optional[VBO] = None
        self._vbo_chords:   Optional[VBO] = None
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
        for attr in ("_vbo_helix", "_vbo_axis_pts", "_vbo_axis_ln",
                     "_vbo_strands", "_vbo_spokes", "_vbo_chords"):
            v = getattr(self, attr)
            if v is not None:
                v.delete()
                setattr(self, attr, None)

    def _rebuild(self) -> None:
        self._free_vbos()
        sps  = generate_semiprimes(self._count)
        span = self._radius * 4.0
        assign_positions(sps, self._radius, span)
        self._stats = compute_stats(sps)

        y_vals = [s.y for s in sps]
        self._y_center = (max(y_vals) + min(y_vals)) / 2 if y_vals else 0.0

        pid_pt = self._pt_prog.id
        pid_ln = self._ln_prog.id

        # Points (stride=6: xyz rgb)
        h_data = _upload_points(sps, filt_lopsided=True)
        a_data = _upload_points(sps, filt_lopsided=False)
        self._vbo_helix    = VBO(h_data, stride=6)
        self._vbo_axis_pts = VBO(a_data, stride=6)
        if h_data: self._vbo_helix.bind_point_attrs(pid_pt)
        if a_data: self._vbo_axis_pts.bind_point_attrs(pid_pt)

        # Lines (stride=7: xyz rgba)
        al_data = build_axis_line(span)
        s_data  = build_strands(sps)
        sp_data = build_spokes(sps)
        c_data  = build_chords(sps)
        self._vbo_axis_ln = VBO(al_data, stride=7)
        self._vbo_strands = VBO(s_data,  stride=7)
        self._vbo_spokes  = VBO(sp_data, stride=7)
        self._vbo_chords  = VBO(c_data,  stride=7)
        for vbo in (self._vbo_axis_ln, self._vbo_strands,
                    self._vbo_spokes, self._vbo_chords):
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
        mvp      = _mvp(self._rot_x, self._rot_y, self._zoom, self._y_center, w, h)
        mvp_flat = (gl.GLfloat * 16)(*mvp)

        # Lines ───────────────────────────────────────────────────────────────
        self._ln_prog.use()
        mvp_loc = gl.glGetUniformLocation(self._ln_prog.id, b"mvp")
        gl.glUniformMatrix4fv(mvp_loc, 1, False, mvp_flat)

        # Axis line is always visible
        if self._vbo_axis_ln and self._vbo_axis_ln.count:
            gl.glBindVertexArray(self._vbo_axis_ln.vao)
            gl.glDrawArrays(gl.GL_LINES, 0, self._vbo_axis_ln.count)
            gl.glBindVertexArray(0)

        for vbo, visible in (
            (self._vbo_strands, self._show_strands),
            (self._vbo_spokes,  self._show_spokes),
            (self._vbo_chords,  self._show_chords),
        ):
            if visible and vbo and vbo.count:
                gl.glBindVertexArray(vbo.vao)
                gl.glDrawArrays(gl.GL_LINES, 0, vbo.count)
                gl.glBindVertexArray(0)

        self._ln_prog.stop()

        # Points ──────────────────────────────────────────────────────────────
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
        if self._vbo_axis_pts and self._vbo_axis_pts.count:
            gl.glBindVertexArray(self._vbo_axis_pts.vao)
            gl.glDrawArrays(gl.GL_POINTS, 0, self._vbo_axis_pts.count)
            gl.glBindVertexArray(0)

        self._pt_prog.stop()

        # 2-D overlay ─────────────────────────────────────────────────────────
        gl.glDisable(gl.GL_DEPTH_TEST)
        self._update_labels()
        self._stat_lbl.draw()
        self._help_lbl.draw()

    def _update_labels(self) -> None:
        s   = self._stats
        fam = s.get("fam", {})
        e   = (
            f"[1]strands={'on' if self._show_strands else 'off'}  "
            f"[2]spokes={'on'  if self._show_spokes  else 'off'}  "
            f"[3]chords={'on'  if self._show_chords  else 'off'}"
        )
        self._stat_lbl.text = (
            f"N={s.get('n', 0):,}  lopsided {s.get('lop_pct', 0):.1f}%  "
            f"balanced {s.get('bal_pct', 0):.1f}%  zoom={abs(self._zoom):.0f}  │  "
            f"1x1={fam.get('1x1', 0):,}  1x3={fam.get('1x3', 0):,}  "
            f"3x3={fam.get('3x3', 0):,}  even={fam.get('even', 0):,}  │  {e}"
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
        # Proportional step: feels smooth at any zoom level.
        # sy > 0 → scroll up → move camera forward (zoom in, less negative).
        step = abs(self._zoom) * 0.08 * sy
        self._zoom = min(-10.0, max(-8_000.0, self._zoom + step))

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
            self._count = max(50, self._count - 100)
            self._rebuild()
        elif symbol == K.BRACKETRIGHT:
            cap = len(_ALL_SP_DATA) if _ALL_SP_DATA else 210_000
            self._count = min(cap, self._count + 100)
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
        self._zoom  = self._zoom_default

# ── CLI / main ────────────────────────────────────────────────────────────────

def _parse() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="PrimeHelix 3D lattice visualizer")
    p.add_argument("--n",      type=int,   default=2000,
                   help="Semiprimes to display (sampled from 1M range, ~73%% lopsided)")
    p.add_argument("--radius", type=float, default=80.0, help="Helix radius")
    p.add_argument("--speed",  type=float, default=0.3,  help="Auto-rotation speed")
    return p.parse_args()


def main() -> None:
    args = _parse()
    sps  = generate_semiprimes(args.n)   # builds sieve cache; prints progress
    lop  = sum(1 for s in sps if s.lopsided)
    print(f"Displaying {len(sps):,} semiprimes — {lop / len(sps) * 100:.1f}% lopsided")
    win = LatticeWindow(args)
    pyglet.app.run()


if __name__ == "__main__":
    main()
