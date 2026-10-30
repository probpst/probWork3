"""
deform.py -- prescribe one cardiac cycle of wall motion on the merged blood-pool STL.

The deformation is applied to the VERTICES of a single base mesh, so every frame
has identical connectivity, identical vertex ordering and identical triangle
ordering.  A displacement field can therefore be obtained by simple subtraction
of two frames, which is what a moving-mesh / ALE solver needs.

Motion model (all fields are smooth in space, so the mesh never tears):
  * left ventricle   longitudinal shortening toward a fixed apex
                     + radial contraction + apex-to-base torsion
  * right ventricle  same construction in its own frame
  * left / right atrium   isotropic reservoir-booster volume change
  * aorta / pulmonary artery   radial distension pulse
  * pulmonary veins, cavae, arch branches   essentially rigid

Usage:  python3 deform.py <base.stl> <output_dir> [n_frames] [bpm]
"""
import os, sys, json
import numpy as np
import trimesh
from scipy.interpolate import PchipInterpolator
import params as P
import geom as G


# ----------------------------------------------------------------------
# volume-versus-time curves (fractions of the cycle, t/T = 0 is end-diastole)
# ----------------------------------------------------------------------
def periodic_curve(tt, vv):
    """Shape-preserving periodic interpolant."""
    tt = np.asarray(tt, float); vv = np.asarray(vv, float)
    t3 = np.r_[tt - 1.0, tt, tt + 1.0]
    v3 = np.r_[vv, vv, vv]
    f = PchipInterpolator(t3, v3)
    return lambda s: f(np.mod(np.asarray(s, float), 1.0))


# left ventricle : EDV 133 mL -> ESV 50 mL  (EF 62 %)
LV_V = periodic_curve(
    [0.00, 0.05, 0.10, 0.16, 0.24, 0.31, 0.36, 0.42, 0.47, 0.53, 0.60, 0.72, 0.85, 0.93],
    [133,  133,  118,   92,   66,   53,   50,   50,   72,  100,  112,  119,  122,  130])
# right ventricle : EDV 140 mL -> ESV 62 mL (EF 56 %)
RV_V = periodic_curve(
    [0.00, 0.05, 0.10, 0.16, 0.24, 0.31, 0.36, 0.43, 0.48, 0.54, 0.62, 0.73, 0.85, 0.93],
    [140,  140,  126,  101,   77,   66,   62,   63,   86,  112,  122,  128,  131,  137])
# left atrium : reservoir max at end LV-systole, min after the atrial kick
LA_V = periodic_curve(
    [0.00, 0.08, 0.18, 0.28, 0.37, 0.43, 0.50, 0.58, 0.70, 0.82, 0.90, 0.96],
    [26,   30,   36,   42,   46,   45,   36,   31,   30,   33,   30,   26])
RA_V = periodic_curve(
    [0.00, 0.08, 0.18, 0.28, 0.37, 0.44, 0.51, 0.59, 0.70, 0.82, 0.90, 0.96],
    [29,   33,   40,   47,   52,   51,   41,   35,   34,   37,   34,   29])
# great-artery diameter pulse (ratio to end-diastolic diameter)
AO_D = periodic_curve(
    [0.00, 0.06, 0.12, 0.20, 0.30, 0.38, 0.50, 0.70, 0.90],
    [1.000, 1.004, 1.035, 1.062, 1.048, 1.026, 1.014, 1.005, 1.001])
PA_D = periodic_curve(
    [0.00, 0.06, 0.14, 0.22, 0.32, 0.40, 0.55, 0.75, 0.90],
    [1.000, 1.004, 1.030, 1.050, 1.038, 1.020, 1.010, 1.004, 1.001])

PHASES = [(0.00, 0.05, 'isovolumic contraction'),
          (0.05, 0.36, 'ventricular ejection'),
          (0.36, 0.43, 'isovolumic relaxation'),
          (0.43, 0.58, 'rapid ventricular filling'),
          (0.58, 0.85, 'diastasis'),
          (0.85, 1.00, 'atrial contraction')]


def phase_name(s):
    s = s % 1.0
    for a, b, n in PHASES:
        if a <= s < b:
            return n
    return PHASES[0][2]


# ----------------------------------------------------------------------
def sstep(x):
    x = np.clip(x, 0.0, 1.0)
    return x * x * (3.0 - 2.0 * x)


def prolate_eta(pts, centre, axis, R, C):
    v = pts - centre
    u = v @ axis
    perp = v - np.outer(u, axis)
    r = np.linalg.norm(perp, axis=1)
    return np.sqrt((r / R) ** 2 + (u / C) ** 2), u, perp


def ellipsoid_eta(pts, centre, radii):
    v = (pts - centre) / radii
    return np.linalg.norm(v, axis=1)


def polyline_radial(pts, ctrl, radii, chunk=20000):
    """Distance to a centreline and the outward radial vector, plus the local radius."""
    c = G.resample_spline(ctrl, 400)
    r = G.resample_scalar(ctrl, radii, 400)
    n = len(pts)
    dist = np.empty(n); rad = np.empty(n); vec = np.empty((n, 3))
    for s in range(0, n, chunk):
        q = pts[s:s + chunk]
        d = q[:, None, :] - c[None, :, :]
        dd = np.linalg.norm(d, axis=2)
        j = np.argmin(dd, axis=1)
        idx = np.arange(len(q))
        dist[s:s + chunk] = dd[idx, j]
        rad[s:s + chunk] = r[j]
        vec[s:s + chunk] = d[idx, j]
    return dist, rad, vec


# ----------------------------------------------------------------------
class Deformer:
    def __init__(self, verts, faces):
        self.V0 = np.asarray(verts, float).copy()
        self.F = np.asarray(faces)
        V = self.V0
        d = P.D_LV

        eta, u, perp = prolate_eta(V, P.O_LV, d, P.LV_R, P.LV_C)
        self.w_lv = sstep((1.44 - eta) / 0.42)
        self.u_lv, self.perp_lv = u, perp

        eta_r, u_r, perp_r = prolate_eta(V, P.O_RV, d, P.RV_R, P.RV_C)
        self.w_rv = sstep((1.32 - eta_r) / 0.30) * (1.0 - self.w_lv)
        self.u_rv, self.perp_rv = u_r, perp_r

        self.w_la = sstep((1.34 - ellipsoid_eta(V, P.LA_C, P.LA_R)) / 0.34)
        self.w_ra = sstep((1.34 - ellipsoid_eta(V, P.RA_C, P.RA_R)) / 0.34)

        dist, rad, vec = polyline_radial(V, P.AORTA_PTS, P.AORTA_R)
        self.w_ao = sstep((1.55 - dist / np.maximum(rad, 1e-6)) / 0.55)
        self.vec_ao = vec
        d1, r1, v1 = polyline_radial(V, P.PA_PTS, P.PA_R)
        w1 = sstep((1.55 - d1 / np.maximum(r1, 1e-6)) / 0.55)
        d2, r2, v2 = polyline_radial(V, P.RPA_PTS, P.RPA_R)
        w2 = sstep((1.55 - d2 / np.maximum(r2, 1e-6)) / 0.55)
        d3, r3, v3 = polyline_radial(V, P.LPA_PTS, P.LPA_R)
        w3 = sstep((1.55 - d3 / np.maximum(r3, 1e-6)) / 0.55)
        self.w_pa = np.maximum(np.maximum(w1, w2), w3)
        sel = np.argmax(np.stack([w1, w2, w3]), axis=0)
        self.vec_pa = np.where(sel[:, None] == 0, v1, np.where(sel[:, None] == 1, v2, v3))

        ss = np.linspace(0, 1, 500)
        self.LV0, self.RV0 = float(LV_V(0.0)), float(RV_V(0.0))
        self.LAmax, self.RAmax = float(LA_V(ss).max()), float(RA_V(ss).max())
        self.a_lv = self.a_rv = 1.0

        # the two disconnected solids (left / right heart); connectivity is fixed
        m = trimesh.Trimesh(V, self.F, process=False)
        comp = trimesh.graph.connected_components(
            m.face_adjacency, nodes=np.arange(len(self.F)))
        probe = np.array([P.O_LV, P.RA_C])          # inside LV / inside RA
        il = 0
        for q, c in enumerate(comp):
            sub = trimesh.Trimesh(V, self.F[c], process=True)
            if bool(sub.contains(probe[:1])[0]):
                il = q
                break
        self.f_left = self.F[comp[il]]
        self.f_right = self.F[comp[1 - il]]

    @staticmethod
    def vol(V, faces):
        a, b, c = V[faces[:, 0]], V[faces[:, 1]], V[faces[:, 2]]
        return abs(float(np.einsum('ij,ij->i', np.cross(a, b), c).sum() / 6.0)) / 1000.0

    @staticmethod
    def _vent_factors(ratio, k_long=0.24):
        ratio = max(ratio, 0.12)
        fL = 1.0 - k_long * (1.0 - ratio)
        return fL, float(np.sqrt(max(ratio / fL, 1e-6)))

    def frame(self, s, parts=('lv', 'rv', 'atria', 'vessels')):
        V = self.V0
        D = np.zeros_like(V)
        d = P.D_LV

        if 'lv' in parts:
            ratio = 1.0 - self.a_lv * (1.0 - float(LV_V(s)) / self.LV0)
            fL, fR = self._vent_factors(ratio)
            zeta = P.LV_C - self.u_lv
            amp = np.clip((1.0 - ratio) / 0.624, 0.0, 1.6)
            phi = np.radians(6.0) * amp * (1.0 - 2.0 * zeta / (P.LV_C + P.LV_BASE_OFF))
            e1 = np.cross(d, np.array([0.0, 0.0, 1.0])); e1 /= np.linalg.norm(e1)
            e2 = np.cross(d, e1)
            a = self.perp_lv @ e1
            b = self.perp_lv @ e2
            ca, sa = np.cos(phi), np.sin(phi)
            a2, b2 = fR * (a * ca - b * sa), fR * (a * sa + b * ca)
            newp = a2[:, None] * e1[None, :] + b2[:, None] * e2[None, :]
            D += self.w_lv[:, None] * ((1.0 - fL) * zeta[:, None] * d[None, :]
                                       + (newp - self.perp_lv))
        if 'rv' in parts:
            ratio = 1.0 - self.a_rv * (1.0 - float(RV_V(s)) / self.RV0)
            fL, fR = self._vent_factors(ratio, 0.20)
            zeta = P.RV_C - self.u_rv
            D += self.w_rv[:, None] * ((1.0 - fL) * zeta[:, None] * d[None, :]
                                       + (fR - 1.0) * self.perp_rv)
        if 'atria' in parts:
            sLA = (float(LA_V(s)) / self.LAmax) ** (1.0 / 3.0)
            D += self.w_la[:, None] * (sLA - 1.0) * (V - P.LA_C[None, :])
            sRA = (float(RA_V(s)) / self.RAmax) ** (1.0 / 3.0)
            D += self.w_ra[:, None] * (sRA - 1.0) * (V - P.RA_C[None, :])
        if 'vessels' in parts:
            D += self.w_ao[:, None] * (float(AO_D(s)) - 1.0) * self.vec_ao
            D += self.w_pa[:, None] * (float(PA_D(s)) - 1.0) * self.vec_pa
        return V + D

    def calibrate(self):
        ss = np.linspace(0, 1, 500)
        for tag, curve, V0c, fset in (('lv', LV_V, self.LV0, self.f_left),
                                      ('rv', RV_V, self.RV0, self.f_right)):
            s_es = float(ss[np.argmin(curve(ss))])
            tgt = float(curve(s_es)) - V0c
            Vrest = self.vol(self.V0, fset)
            lo, hi = 0.6, 1.8
            for _ in range(26):
                mid = 0.5 * (lo + hi)
                setattr(self, 'a_' + tag, mid)
                dv = self.vol(self.frame(s_es, parts=(tag,)), fset) - Vrest
                if dv > tgt:
                    lo = mid
                else:
                    hi = mid
            setattr(self, 'a_' + tag, 0.5 * (lo + hi))
            dv = self.vol(self.frame(s_es, parts=(tag,)), fset) - Vrest
            print(f"  calibrate {tag}: amplitude x{getattr(self,'a_'+tag):.3f}  "
                  f"target {tgt:+.1f} mL  achieved {dv:+.1f} mL")


# ----------------------------------------------------------------------
if __name__ == '__main__':
    base = sys.argv[1] if len(sys.argv) > 1 else '/home/claude/cycle_base.stl'
    out = sys.argv[2] if len(sys.argv) > 2 else '/home/claude/cycle_out'
    nfr = int(sys.argv[3]) if len(sys.argv) > 3 else 40
    bpm = float(sys.argv[4]) if len(sys.argv) > 4 else 70.0
    T = 60.0 / bpm

    os.makedirs(out, exist_ok=True)
    M = trimesh.load(base)
    dfm = Deformer(M.vertices, M.faces)
    dfm.calibrate()
    meta = dict(base_stl=os.path.basename(base), n_frames=nfr, bpm=bpm,
                cycle_period_s=round(T, 5), units='mm',
                n_triangles=int(len(M.faces)), n_vertices=int(len(M.vertices)),
                note='identical connectivity and vertex ordering in every frame',
                frames=[])
    print(f"{'frame':>5} {'t [s]':>8} {'t/T':>6} {'phase':26s} "
          f"{'V_left':>8} {'V_right':>8} {'minTriA':>9} {'maxDisp':>8}")
    prev = None
    for n in range(nfr):
        s = n / nfr
        Vn = dfm.frame(s)
        m = trimesh.Trimesh(Vn, M.faces, process=False)
        vols = [dfm.vol(Vn, dfm.f_left), dfm.vol(Vn, dfm.f_right)]
        tri = m.area_faces
        disp = 0.0 if prev is None else float(np.linalg.norm(Vn - prev, axis=1).max())
        fn = os.path.join(out, f"heart_cycle_{n:03d}.stl")
        m.export(fn)
        print(f"{n:5d} {s*T:8.4f} {s:6.3f} {phase_name(s):26s} "
              f"{vols[0]:8.1f} {vols[1]:8.1f} "
              f"{tri.min():9.2e} {disp:8.3f}")
        meta['frames'].append(dict(frame=n, file=os.path.basename(fn),
                                   t_s=round(s * T, 5), t_over_T=round(s, 4),
                                   phase=phase_name(s),
                                   LV_mL=round(float(LV_V(s)), 1),
                                   RV_mL=round(float(RV_V(s)), 1),
                                   LA_mL=round(float(LA_V(s)), 1),
                                   RA_mL=round(float(RA_V(s)), 1),
                                   aortic_diam_ratio=round(float(AO_D(s)), 4),
                                   pulm_diam_ratio=round(float(PA_D(s)), 4),
                                   total_volume_mL=round(float(sum(vols)), 2),
                                   max_vertex_step_mm=round(disp, 4)))
        prev = Vn
    with open(os.path.join(out, 'cycle_metadata.json'), 'w') as f:
        json.dump(meta, f, indent=2)
    print("\nwritten ->", out)
