"""Fig. 1: model cross-section and winding layout (run from analysis/)."""
import sys, math
sys.path.insert(0, '../fea')
import im_geometry as geo
geo.R_ROTOR = 62.8 - 0.44
from shapely.geometry import LineString, Point
from shapely.ops import unary_union, polygonize
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch, Patch
nodes, segs, arcs, labels = geo.build()
rd = lambda p: (round(p[0], 6), round(p[1], 6))
L = [LineString([rd(p), rd(q)]) for p, q in segs]
for p, q, sw in arcs:
    a = math.radians(sw); ch = math.dist(p, q); r = ch / (2 * math.sin(a / 2)); mx, my = (p[0] + q[0]) / 2, (p[1] + q[1]) / 2
    dx, dy = q[0] - p[0], q[1] - p[1]; h = math.sqrt(max(r * r - ch * ch / 4, 0)) * (1 if sw <= 180 else -1)
    cx, cy = mx - dy / ch * h, my + dx / ch * h; a0 = math.atan2(p[1] - cy, p[0] - cx); n = max(6, int(sw / 0.5))
    pts = [(cx + r * math.cos(a0 + a * i / n), cy + r * math.sin(a0 + a * i / n)) for i in range(n + 1)]; pts[0] = rd(p); pts[-1] = rd(q)
    L.append(LineString(pts))
F = list(polygonize(unary_union(L)))
BELTS = [("A", 1), ("C", -1), ("B", 1), ("A", -1), ("C", 1), ("B", -1)]
up = lambda k: BELTS[(k // 3) % 6]; lo = lambda k: (up((k - 7) % 36)[0], -up((k - 7) % 36)[1])
PC = {'A': '#2a78d6', 'B': '#eb6834', 'C': '#1baf7a'}
def path_of(poly):
    verts, codes = [], []
    for ring in [poly.exterior] + list(poly.interiors):
        c = list(ring.coords); verts += c; codes += [Path.MOVETO] + [Path.LINETO] * (len(c) - 2) + [Path.CLOSEPOLY]
    return Path(verts, codes)
fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.6), gridspec_kw={'width_ratios': [1, 1.05]})
for ax, zoom in zip(axs, (False, True)):
    for f in F:
        col, hatch = '#ffffff', None
        for role, x, y, k in labels:
            if f.contains(Point(x, y)):
                if role in ('stator_iron', 'rotor_iron'): col = '#d9d8d2'
                elif role == 'shaft': col = '#b8b7ae'
                elif role == 'bar': col = '#8f8f8a'
                elif role in ('coil_up', 'coil_lo'):
                    ph, sg = up(k) if role == 'coil_up' else lo(k); col = PC[ph]; hatch = None if sg > 0 else '////'
                break
        ax.add_patch(PathPatch(path_of(f), fc=col, ec='#22221f', lw=0.15 if not zoom else 0.5, hatch=hatch))
    ax.set_aspect('equal'); ax.axis('off')
    if zoom: ax.set_xlim(38, 82); ax.set_ylim(-2, 28)
    else: ax.set_xlim(-101, 101); ax.set_ylim(-101, 101)
axs[0].set_title('(a) 2D model, 36 slots / 28 bars', fontsize=9, loc='left', family='DejaVu Serif')
axs[1].set_title('(b) Slot detail: double layer, pitch 7/9', fontsize=9, loc='left', family='DejaVu Serif')
h = [Patch(fc=PC[p], ec='#22221f', label='phase %s' % {'A': 'U', 'B': 'V', 'C': 'W'}[p]) for p in 'ABC'] + \
    [Patch(fc='white', ec='#22221f', hatch='////', label='return side (−)'), Patch(fc='#8f8f8a', ec='#22221f', label='Al bar')]
axs[1].legend(handles=h, fontsize=7, loc='upper center', bbox_to_anchor=(0.5, -0.02), ncol=3, frameon=False)
fig.savefig('../figures/fig_model.png', dpi=300, bbox_inches='tight'); fig.savefig('../figures/fig_model.pdf', bbox_inches='tight')
