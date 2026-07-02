#!/usr/bin/env python3
"""Generate Fig. 6 (Poinsot geometry) for index.html.

Every curve is computed from the actual quadric equations and projected
orthographically, so the drawing cannot disagree with the physics. Each
sampled point is asserted to lie on both quadrics before it is drawn.

Usage:  python tools/gen_fig6.py
Replaces the markup between the FIG6-GENERATED markers in index.html.

ponytail: silhouettes are drawn without mutual occlusion between the two
bodies in the left panel (textbook style); add ray-tested clipping if a
camera angle ever makes it misleading.
"""
import math
import os

# principal moments, shared by both panels (same values as Exercise B.3)
I1, I2, I3 = 1.0, 2.0, 3.0

# ---------------------------------------------------------------- vectors
def dot(a, b):
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]

def cross(a, b):
    return (a[1]*b[2] - a[2]*b[1],
            a[2]*b[0] - a[0]*b[2],
            a[0]*b[1] - a[1]*b[0])

# ---------------------------------------------------------------- camera
AZ, EL = math.radians(35.0), math.radians(18.0)
VIEW  = (math.cos(EL)*math.cos(AZ), math.cos(EL)*math.sin(AZ), math.sin(EL))
RIGHT = (-math.sin(AZ), math.cos(AZ), 0.0)
UP    = cross(VIEW, RIGHT)          # unit because VIEW ⟂ RIGHT

def project(p, scale, cx, cy):
    return (cx + scale*dot(p, RIGHT), cy - scale*dot(p, UP))

def fmt(xy):
    return "%.1f %.1f" % xy

def path_d(pts, close=False):
    d = "M" + " L".join(fmt(p) for p in pts)
    return d + " Z" if close else d

def visibility_paths(pts3, visible, scale, cx, cy):
    """Split a closed 3D curve into visible/hidden projected polylines."""
    n = len(pts3)
    vis_ds, hid_ds = [], []
    cur, cur_flag = [], None
    for k in range(n + 1):
        i = k % n
        f, p = visible[i], project(pts3[i], scale, cx, cy)
        if cur_flag is None:
            cur, cur_flag = [p], f
        elif f == cur_flag:
            cur.append(p)
        else:
            cur.append(p)                      # share the transition point
            (vis_ds if cur_flag else hid_ds).append(path_d(cur))
            cur, cur_flag = [p], f
    if len(cur) > 1:
        (vis_ds if cur_flag else hid_ds).append(path_d(cur))
    return vis_ds, hid_ds

def silhouette_d(semi, scale, cx, cy, n=180):
    """Orthographic outline of ellipsoid diag(1/a_i^2): y^T (B^T A^{-1} B)^{-1} y = 1."""
    a2 = [s*s for s in semi]
    n11 = sum(a2[i]*RIGHT[i]*RIGHT[i] for i in range(3))
    n12 = sum(a2[i]*RIGHT[i]*UP[i]    for i in range(3))
    n22 = sum(a2[i]*UP[i]*UP[i]       for i in range(3))
    tr, det = n11 + n22, n11*n22 - n12*n12
    disc = math.sqrt(max(tr*tr/4.0 - det, 0.0))
    l1, l2 = tr/2.0 + disc, tr/2.0 - disc
    if abs(n12) > 1e-12:
        v1 = (l1 - n22, n12)
    else:
        v1 = (1.0, 0.0) if n11 >= n22 else (0.0, 1.0)
    m = math.hypot(*v1)
    v1 = (v1[0]/m, v1[1]/m)
    v2 = (-v1[1], v1[0])
    r1, r2 = math.sqrt(l1), math.sqrt(l2)
    pts = []
    for k in range(n):
        t = 2.0*math.pi*k/n
        u = r1*math.cos(t)*v1[0] + r2*math.sin(t)*v2[0]
        w = r1*math.cos(t)*v1[1] + r2*math.sin(t)*v2[1]
        pts.append((cx + scale*u, cy - scale*w))
    return path_d(pts, close=True)

# ============================== LEFT PANEL ==============================
# Momentum space: sphere |Pi|^2 = L2 intersected with energy ellipsoid
# sum Pi_i^2/(2 I_i) = T.  With T < L2/(2 I2) the curve encircles axis 3.
L2, T = 1.0, 0.22
A = (1.0/(2*I1), 1.0/(2*I2), 1.0/(2*I3))
assert A[2]*L2 < T < A[1]*L2, "energy must put the curve around axis 3"
P_AMP = math.sqrt((T - A[2]*L2)/(A[0] - A[2]))
Q_AMP = math.sqrt((T - A[2]*L2)/(A[1] - A[2]))
E_SEMI = tuple(math.sqrt(2*T*Ii) for Ii in (I1, I2, I3))

def pi_curve(sign, n=240):
    pts = []
    for k in range(n):
        t = 2.0*math.pi*k/n
        x, y = P_AMP*math.cos(t), Q_AMP*math.sin(t)
        z2 = L2 - x*x - y*y
        assert z2 > 0.0
        pt = (x, y, sign*math.sqrt(z2))
        assert abs(dot(pt, pt) - L2) < 1e-9
        assert abs(sum(Ai*c*c for Ai, c in zip(A, pt)) - T) < 1e-9
        pts.append(pt)
    return pts

def pi_visible(pt):
    en = (pt[0]/E_SEMI[0]**2, pt[1]/E_SEMI[1]**2, pt[2]/E_SEMI[2]**2)
    return dot(pt, VIEW) > 0.0 and dot(en, VIEW) > 0.0

# ============================== RIGHT PANEL =============================
# Omega space with 2T = 1: inertia ellipsoid sum I_i w_i^2 = 1; polhodes are
# its intersections with sum I_i^2 w_i^2 = L2p for L2p in (I1, I3).
E = (I1, I2, I3)
B_SEMI = tuple(1.0/math.sqrt(Ii) for Ii in E)

def polhode(k, L2p, sign, n=240):
    i, j = [m for m in range(3) if m != k]
    ci, cj = E[i]*(E[i] - E[k]), E[j]*(E[j] - E[k])
    rhs = L2p - E[k]
    ri, rj = math.sqrt(rhs/ci), math.sqrt(rhs/cj)
    pts = []
    for s in range(n):
        t = 2.0*math.pi*s/n
        w = [0.0, 0.0, 0.0]
        w[i], w[j] = ri*math.cos(t), rj*math.sin(t)
        rem = 1.0 - E[i]*w[i]**2 - E[j]*w[j]**2
        assert rem > 0.0
        w[k] = sign*math.sqrt(rem/E[k])
        pt = tuple(w)
        assert abs(sum(Ei*c*c for Ei, c in zip(E, pt)) - 1.0) < 1e-9
        assert abs(sum(Ei*Ei*c*c for Ei, c in zip(E, pt)) - L2p) < 1e-9
        pts.append(pt)
    return pts

def separatrix(sign, n=240):
    slope = sign*math.sqrt((E[0]*(E[1] - E[0]))/(E[2]*(E[2] - E[1])))
    r1 = 1.0/math.sqrt(E[0] + E[2]*slope*slope)
    r2 = 1.0/math.sqrt(E[1])
    pts = []
    for s in range(n):
        t = 2.0*math.pi*s/n
        x1 = r1*math.cos(t)
        pt = (x1, r2*math.sin(t), slope*x1)
        assert abs(sum(Ei*c*c for Ei, c in zip(E, pt)) - 1.0) < 1e-9
        assert abs(sum(Ei*Ei*c*c for Ei, c in zip(E, pt)) - E[1]) < 1e-9
        pts.append(pt)
    return pts

def omega_visible(pt):
    return dot((E[0]*pt[0], E[1]*pt[1], E[2]*pt[2]), VIEW) > 0.0

# ============================== SVG ASSEMBLY ============================
LCX, LCY, LS = 180, 160, 105     # left panel center and scale
RCX, RCY, RS = 515, 160, 140     # right panel center and scale

def curve_svg(pts3, visible_fn, scale, cx, cy, color, width):
    vis = [visible_fn(p) for p in pts3]
    vis_ds, hid_ds = visibility_paths(pts3, vis, scale, cx, cy)
    out = []
    for d in hid_ds:
        out.append('<path d="%s" fill="none" stroke="%s" stroke-width="%.1f" '
                   'stroke-dasharray="5 5" opacity=".45"></path>' % (d, color, width*0.75))
    for d in vis_ds:
        out.append('<path d="%s" fill="none" stroke="%s" stroke-width="%.1f"></path>'
                   % (d, color, width))
    return out

parts = []
parts.append('<svg viewBox="0 0 700 320" width="100%" role="img" '
             'aria-label="exact Poinsot construction: momentum sphere and energy '
             'ellipsoid intersection, and polhode family on the inertia ellipsoid">')
parts.append('<rect x="0" y="0" width="700" height="320" fill="#fff"></rect>')

# ---- left: sphere silhouette, ellipsoid silhouette, intersection curve
parts.append('<circle cx="%d" cy="%d" r="%.1f" fill="#f8f5ec" stroke="#1c5f74" '
             'stroke-width="2.2"></circle>' % (LCX, LCY, LS*math.sqrt(L2)))
parts.append('<path d="%s" fill="none" stroke="#c08a2d" stroke-width="2"></path>'
             % silhouette_d(E_SEMI, LS, LCX, LCY))
top = pi_curve(+1)
parts += curve_svg(top, pi_visible, LS, LCX, LCY, "#7a1f1f", 2.8)
parts += curve_svg(pi_curve(-1), pi_visible, LS, LCX, LCY, "#7a1f1f", 2.8)
# animated dot along the (closed) upper branch
parts.append('<circle r="5" fill="#7a1f1f"><animateMotion dur="8s" '
             'repeatCount="indefinite" path="%s"></animateMotion></circle>'
             % path_d([project(p, LS, LCX, LCY) for p in top], close=True))
# label the curve at its rightmost projected point
lx, ly = max((project(p, LS, LCX, LCY) for p in top), key=lambda q: q[0])
parts.append('<text class="v" x="%.0f" y="%.0f" font-size="13" fill="#7a1f1f">'
             'path of &#928;</text>' % (lx + 8, ly - 4))
parts.append('<text x="52" y="34" font-size="13" fill="#1c5f74">momentum sphere</text>')
parts.append('<text x="60" y="306" font-size="13" fill="#c08a2d">energy ellipsoid</text>')

# ---- right: inertia ellipsoid silhouette, axes, polhodes, separatrix
parts.append('<path d="%s" fill="#f6f4ee" stroke="#9ca9b1" stroke-width="2.2"></path>'
             % silhouette_d(B_SEMI, RS, RCX, RCY))
for k, sub in ((0, "&#8321;"), (1, "&#8322;"), (2, "&#8323;")):
    tip = [0.0, 0.0, 0.0]
    tip[k] = 1.28*B_SEMI[k]
    x0, y0 = project(tuple(-t for t in tip), RS, RCX, RCY)
    x1, y1 = project(tuple(tip), RS, RCX, RCY)
    parts.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#999" '
                 'stroke-width="1.4" stroke-dasharray="4 4"></line>' % (x0, y0, x1, y1))
    parts.append('<text class="v" x="%.0f" y="%.0f" font-size="13" fill="#555">'
                 '&#969;<tspan baseline-shift="sub" font-size="10">%s</tspan></text>'
                 % (x1 + 5, y1 - 3, sub))
for L2p in (1.12, 1.5):                       # around stable axis 1
    for sign in (+1, -1):
        parts += curve_svg(polhode(0, L2p, sign), omega_visible, RS, RCX, RCY, "#7a1f1f", 2.2)
for L2p in (2.5, 2.88):                       # around stable axis 3
    for sign in (+1, -1):
        parts += curve_svg(polhode(2, L2p, sign), omega_visible, RS, RCX, RCY, "#7a1f1f", 2.2)
for sign in (+1, -1):                         # separatrix through axis 2
    parts += curve_svg(separatrix(sign), omega_visible, RS, RCX, RCY, "#2f6d45", 2.4)
parts.append('<text x="404" y="34" font-size="13" fill="#7a1f1f">polhodes on inertia ellipsoid</text>')
parts.append('<text x="480" y="306" font-size="13" fill="#2f6d45">separatrix</text>')
parts.append('</svg>')
svg = "\n".join(parts)

# ============================== SPLICE INTO HTML ========================
BEGIN, END = "<!-- FIG6-GENERATED:BEGIN", "<!-- FIG6-GENERATED:END -->"
html_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "index.html")
with open(html_path, encoding="utf-8") as f:
    html = f.read()
b = html.index(BEGIN)
b_end = html.index("-->", b) + 3
e = html.index(END)
assert html.count(BEGIN) == 1 and html.count(END) == 1
html = html[:b_end] + "\n" + svg + "\n" + html[e:]
with open(html_path, "w", encoding="utf-8", newline="\n") as f:
    f.write(html)
print("Fig. 6 regenerated: %d SVG elements, all curve points verified on both quadrics."
      % svg.count("<path"))
