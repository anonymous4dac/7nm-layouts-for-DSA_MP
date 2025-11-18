import gdspy
import numpy as np
from dataclasses import dataclass
from typing import Tuple, Optional, Dict, List
from scipy.spatial.distance import cdist
import networkx as nx
from collections import Counter
from networkx.algorithms.clique import find_cliques

# ============================
# Metal Layer Specification
# ============================

@dataclass
class MetalGenSpec:
    wire_cd: float
    track_pitch: float
    min_t2t: float
    max_t2t: float
    t2t_grid: float
    min_length: float
    max_length: float
    total_x: float
    total_y: float
    cell_name: str = "metal_cell"
    layer: int = 1
    datatype: int = 0
    origin: Tuple[float, float] = (0.0, 0.0)
    seed: Optional[int] = 0


# ============================
# Helper Functions
# ============================

def _rand_quantized(rng, lo, hi, grid):
    if hi < lo:
        hi, lo = lo, hi
    if grid <= 0:
        return float(rng.uniform(lo, hi))
    val = rng.uniform(lo, hi)
    q = round(val / grid) * grid
    return min(max(q, lo), hi)


def _draw_one_track_vertical(cell, x_base, spec, rng):
    w = spec.wire_cd
    yt = spec.total_y
    y = 0.0
    xo, yo = spec.origin

    while y < yt:
        rem = yt - y
        l_hi = min(spec.max_length, rem)
        if l_hi < spec.min_length:
            break

        l = _rand_quantized(rng, spec.min_length, l_hi, 0.0)

        rect = gdspy.Rectangle(
            (xo + x_base, yo + y),
            (xo + x_base + w, yo + y + l),
            layer=spec.layer,
            datatype=spec.datatype
        )
        cell.add(rect)

        y += l
        if y >= yt:
            break

        t_hi = min(spec.max_t2t, yt - y)
        if t_hi < spec.min_t2t:
            break

        t = _rand_quantized(rng, spec.min_t2t, t_hi, spec.t2t_grid)
        y += t


def draw_metal_cell_vertical(spec: MetalGenSpec) -> gdspy.GdsLibrary:
    lib = gdspy.GdsLibrary(unit=1e-9, precision=1e-12)
    cell = lib.new_cell(spec.cell_name)
    rng = np.random.default_rng(spec.seed)

    x = 0.0
    while x + spec.wire_cd <= spec.total_x:
        _draw_one_track_vertical(cell, x, spec, rng)
        x += spec.track_pitch

    return lib

def draw_metal1(
    lib: gdspy.GdsLibrary,
    length_x: float,
    width: float,
    pitch: float,
    total_y: float,
    layer: int = 110,
    datatype: int = 0,
) -> gdspy.GdsLibrary:
    cell = list(lib.cells.values())[0]
    y = 0.0
    while y + width <= total_y:
        rect = gdspy.Rectangle(
            (0.0, y),
            (length_x, y + width),
            layer=layer,
            datatype=datatype,
        )
        cell.add(rect)
        y += pitch
    return lib


def save_gds(lib, out_path):
    lib.write_gds(out_path)


# ============================
# Via Specification
# ============================

@dataclass
class ViaGenSpec:
    via_x: float
    via_y: float
    density: float
    enclosure_x: float
    enclosure_y: float
    via_pitch_x: float
    via_pitch_y: float
    via_layer: int = 112
    via_datatype: int = 0
    seed: int = 0


# ============================
# GDS Geometry Helpers
# ============================

def _get_layer_polys(cell, layer, datatype=0):
    by_spec = cell.get_polygons(by_spec=True)
    return by_spec.get((layer, datatype), [])


def _bbox(poly):
    xs, ys = poly[:, 0], poly[:, 1]
    return float(xs.min()), float(ys.min()), float(xs.max()), float(ys.max())


def _rect_from_bbox(b):
    x0, y0, x1, y1 = b
    if x1 <= x0 or y1 <= y0:
        return None
    return gdspy.Rectangle((x0, y0), (x1, y1))


def _unique_positions(vals, tol):
    if not vals:
        return []
    vals = sorted(vals)
    uniq = [vals[0]]
    for v in vals[1:]:
        if abs(v - uniq[-1]) > tol:
            uniq.append(v)
    return uniq


# ============================
# Assist Wires
# ============================

def _assist_polys_horizontal(polys, ex):
    out = []
    for p in polys:
        x0, y0, x1, y1 = _bbox(p)
        r = _rect_from_bbox((x0 + ex, y0, x1 - ex, y1))
        if r:
            out.append(r)
    return out


def _assist_polys_vertical(polys, ey):
    out = []
    for p in polys:
        x0, y0, x1, y1 = _bbox(p)
        r = _rect_from_bbox((x0, y0 + ey, x1, y1 - ey))
        if r:
            out.append(r)
    return out


# ============================
# Via Candidate & Pitch Cleaning
# ============================

def _calc_candidates(cell, m1_layer, m2_layer, ex_nm, ey_nm):
    polys_m1 = _get_layer_polys(cell, m1_layer)
    polys_m2 = _get_layer_polys(cell, m2_layer)
    assist_m1 = _assist_polys_horizontal(polys_m1, ex_nm)
    assist_m2 = _assist_polys_vertical(polys_m2, ey_nm)

    if not assist_m1 or not assist_m2:
        return []

    inter = gdspy.boolean(assist_m1, assist_m2, 'and', precision=1e-12, max_points=4094)
    if inter is None:
        return []

    polys = []
    if isinstance(inter, gdspy.Polygon):
        polys = [inter.points]
    else:
        polys = inter.polygons

    bboxes = [_bbox(p) for p in polys]
    return bboxes


def _sample_centers_from_bboxes(bboxes, vx, vy, rng, density):
    centers = []
    for x0, y0, x1, y1 in bboxes:
        if (x1 - x0) >= vx and (y1 - y0) >= vy:
            if rng.random() < density:
                centers.append(((x0 + x1) * 0.5, (y0 + y1) * 0.5))
    return centers


def _pitch_clean(centers, vpx, vpy, polys_m1, polys_m2):
    if not centers:
        return []

    x_tracks = _unique_positions([(_bbox(p)[0] + _bbox(p)[2]) * 0.5 for p in polys_m2], tol=1e-3)
    y_tracks = _unique_positions([(_bbox(p)[1] + _bbox(p)[3]) * 0.5 for p in polys_m1], tol=1e-3)

    def nearest_idx(v, arr):
        if not arr:
            return -1
        idx = min(range(len(arr)), key=lambda i: abs(arr[i] - v))
        return idx

    # Column cleaning (y-direction pitch)
    col_map = {}
    for cx, cy in centers:
        col = nearest_idx(cx, x_tracks)
        col_map.setdefault(col, []).append((cx, cy))

    kept = []
    for pts in col_map.values():
        pts.sort(key=lambda t: t[1])
        last = None
        for pt in pts:
            if last is None or (pt[1] - last[1]) >= vpy:
                kept.append(pt)
                last = pt

    # Row cleaning (x-direction pitch)
    row_map = {}
    for cx, cy in kept:
        row = nearest_idx(cy, y_tracks)
        row_map.setdefault(row, []).append((cx, cy))

    final = []
    for pts in row_map.values():
        pts.sort(key=lambda t: t[0])
        last = None
        for pt in pts:
            if last is None or (pt[0] - last[0]) >= vpx:
                final.append(pt)
                last = pt
    return final


# ============================
# VIA Generation
# ============================

def generate_vias_from_metals(lib,
                              spec: ViaGenSpec,
                              m1_layer=110,
                              m2_layer=111):

    cell = list(lib.cells.values())[0]
    rng = np.random.default_rng(spec.seed)

    bboxes = _calc_candidates(cell, m1_layer, m2_layer,
                              ex_nm=spec.enclosure_x,
                              ey_nm=spec.enclosure_y)

    centers = _sample_centers_from_bboxes(
        bboxes,
        vx=spec.via_x,
        vy=spec.via_y,
        rng=rng,
        density=spec.density
    )

    polys_m1 = _get_layer_polys(cell, m1_layer)
    polys_m2 = _get_layer_polys(cell, m2_layer)
    centers = _pitch_clean(
        centers,
        vpx=spec.via_pitch_x,
        vpy=spec.via_pitch_y,
        polys_m1=polys_m1,
        polys_m2=polys_m2
    )

    vias = []
    hx, hy = 0.5 * spec.via_x, 0.5 * spec.via_y
    for cx, cy in centers:
        rect = gdspy.Rectangle(
            (cx - hx, cy - hy),
            (cx + hx, cy + hy),
            layer=spec.via_layer,
            datatype=spec.via_datatype
        )
        cell.add(rect)
        vias.append(rect)

    return vias


# ============================
# Example Usage
# ============================

if __name__ == "__main__":

    spec = MetalGenSpec(
        wire_cd=21,
        track_pitch=31.5,
        min_t2t=42,
        max_t2t=315,
        t2t_grid=21,
        min_length=42,
        max_length=378,
        total_x=5533,
        total_y=5187,
        layer=111
    )

    lib = draw_metal_cell_vertical(spec)

    lib = draw_metal1(
        lib,
        length_x=spec.total_x,
        width=21,
        pitch=42,
        total_y=spec.total_y,
        layer=110          # 横向 metal 层
    )

    via_spec = ViaGenSpec(
        via_x=21,
        via_y=21,
        density=0.5,
        enclosure_x=5,
        enclosure_y=5,
        via_pitch_x=32.5,
        via_pitch_y=22
    )

    generate_vias_from_metals(lib, via_spec)

    # Generic output path for GitHub projects
    save_gds(lib, "generated_layout.gds")
