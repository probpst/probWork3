"""
geom.py -- parametric surface construction utilities for cardiovascular CFD geometry.
All units in millimetres.
Coordinate system (patient / body frame):
    +x : LEFT   (patient's left)
    +y : ANTERIOR (front of chest)
    +z : SUPERIOR (height above floor)  <-- z is the height axis
    origin : mid-sagittal plane, mid-coronal plane, floor level
"""
import numpy as np
import trimesh
from scipy.interpolate import CubicSpline


# ----------------------------------------------------------------------
# centreline helpers
# ----------------------------------------------------------------------
def resample_spline(ctrl, n=200):
    """Cubic-spline interpolation of control points, resampled uniformly in arclength."""
    ctrl = np.asarray(ctrl, float)
    # chord-length parameterisation
    d = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(ctrl, axis=0), axis=1))]
    t = d / d[-1]
    cs = CubicSpline(t, ctrl, axis=0, bc_type='natural')
    tt = np.linspace(0, 1, 4000)
    dense = cs(tt)
    dd = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(dense, axis=0), axis=1))]
    s = np.linspace(0, dd[-1], n)
    out = np.empty((n, 3))
    for k in range(3):
        out[:, k] = np.interp(s, dd, dense[:, k])
    return out


def resample_scalar(ctrl_pts, ctrl_val, n):
    """Interpolate a scalar (e.g. radius) defined at the control points onto n samples."""
    ctrl_pts = np.asarray(ctrl_pts, float)
    d = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(ctrl_pts, axis=0), axis=1))]
    t = d / d[-1]
    cs = CubicSpline(t, np.asarray(ctrl_val, float), bc_type='natural')
    return cs(np.linspace(0, 1, n))


def parallel_frames(curve):
    """Rotation-minimising (parallel transport) frames along a curve.
    Returns tangent T, normal U, binormal V arrays of shape (n,3)."""
    n = len(curve)
    T = np.gradient(curve, axis=0)
    T /= np.linalg.norm(T, axis=1)[:, None]
    # seed normal
    a = np.array([0.0, 0.0, 1.0])
    if abs(np.dot(a, T[0])) > 0.9:
        a = np.array([1.0, 0.0, 0.0])
    U = np.empty_like(T)
    U[0] = np.cross(T[0], a)
    U[0] /= np.linalg.norm(U[0])
    for i in range(1, n):
        v = np.cross(T[i - 1], T[i])
        s = np.linalg.norm(v)
        if s < 1e-9:
            U[i] = U[i - 1]
        else:
            v /= s
            ang = np.arctan2(s, np.dot(T[i - 1], T[i]))
            c, si = np.cos(ang), np.sin(ang)
            u = U[i - 1]
            U[i] = u * c + np.cross(v, u) * si + v * np.dot(v, u) * (1 - c)
        U[i] -= T[i] * np.dot(U[i], T[i])
        U[i] /= np.linalg.norm(U[i])
    V = np.cross(T, U)
    return T, U, V


# ----------------------------------------------------------------------
# ring / loft primitives
# ----------------------------------------------------------------------
def circle_profile(nphi):
    ang = np.linspace(0, 2 * np.pi, nphi, endpoint=False)
    return np.cos(ang), np.sin(ang)


def sweep_rings(curve, radii, nphi=96, profile=None, twist=None):
    """Build rings of points around a centreline.
    profile: optional callable(i, ang) -> radial scale factor (1.0 = circular)."""
    T, U, V = parallel_frames(curve)
    ang = np.linspace(0, 2 * np.pi, nphi, endpoint=False)
    ca, sa = np.cos(ang), np.sin(ang)
    rings = np.empty((len(curve), nphi, 3))
    for i in range(len(curve)):
        sc = 1.0
        if profile is not None:
            sc = profile(i, ang)
        r = radii[i] * sc
        rings[i] = curve[i][None, :] + (r * ca)[:, None] * U[i][None, :] \
                                     + (r * sa)[:, None] * V[i][None, :]
    return rings


def loft_faces(nsec, nphi, closed_ring=True):
    """Triangle connectivity for a stack of nsec rings of nphi points."""
    f = []
    for i in range(nsec - 1):
        for j in range(nphi):
            jn = (j + 1) % nphi if closed_ring else min(j + 1, nphi - 1)
            if not closed_ring and j == nphi - 1:
                continue
            a = i * nphi + j
            b = i * nphi + jn
            c = (i + 1) * nphi + jn
            d = (i + 1) * nphi + j
            f.append([a, b, c])
            f.append([a, c, d])
    return np.asarray(f, np.int64)


def cap_delaunay(ring3d, target_edge=0.8):
    """Well-shaped planar triangulation of a (convex) boundary ring.
    Returns (vertices3d, faces) with the ring occupying the first len(ring) rows."""
    from scipy.spatial import Delaunay
    centre = ring3d.mean(axis=0)
    Q = ring3d - centre
    u, s, vt = np.linalg.svd(Q, full_matrices=False)
    e1, e2, nrm = vt[0], vt[1], vt[2]
    P2 = np.column_stack([Q @ e1, Q @ e2])
    R = np.linalg.norm(P2, axis=1).max()
    h = target_edge
    pts = [P2]
    ny = int(np.ceil(2 * R / (h * np.sqrt(3) / 2))) + 1
    for iy in range(-ny, ny + 1):
        y = iy * h * np.sqrt(3) / 2
        xoff = (h / 2) if (iy % 2) else 0.0
        nx = int(np.ceil(2 * R / h)) + 1
        xs = np.arange(-nx, nx + 1) * h + xoff
        cand = np.column_stack([xs, np.full_like(xs, y)])
        rr = np.linalg.norm(cand, axis=1)
        # radial boundary of the (convex) ring, by angle
        ang = np.arctan2(cand[:, 1], cand[:, 0])
        ra = np.linalg.norm(P2, axis=1)
        aa = np.arctan2(P2[:, 1], P2[:, 0])
        o = np.argsort(aa)
        rb = np.interp(ang, aa[o], ra[o], period=2 * np.pi)
        keep = rr < rb - 0.62 * h
        if keep.any():
            pts.append(cand[keep])
    allp = np.vstack(pts)
    tri = Delaunay(allp)
    f = tri.simplices
    cen = allp[f].mean(axis=1)
    rc = np.linalg.norm(cen, axis=1)
    ac = np.arctan2(cen[:, 1], cen[:, 0])
    ra = np.linalg.norm(P2, axis=1); aa = np.arctan2(P2[:, 1], P2[:, 0])
    o = np.argsort(aa)
    f = f[rc < np.interp(ac, aa[o], ra[o], period=2 * np.pi) + 1e-9]
    V3 = centre + allp[:, 0:1] * e1[None, :] + allp[:, 1:2] * e2[None, :]
    return V3, f


def cap_fan(ring_pts, base_index, reverse=False):
    """Fan-triangulate a ring using its centroid. Returns (extra_vertex, faces)."""
    nphi = len(ring_pts)
    centre = ring_pts.mean(axis=0)
    ci = base_index + nphi  # centroid appended after the ring
    f = []
    for j in range(nphi):
        jn = (j + 1) % nphi
        if reverse:
            f.append([base_index + jn, base_index + j, ci])
        else:
            f.append([base_index + j, base_index + jn, ci])
    return centre, np.asarray(f, np.int64)


def _attach_cap(rings_edge, base_index, verts, faces, nv, outward, target_edge):
    """Append a Delaunay-triangulated cap on the given boundary ring."""
    nphi = len(rings_edge)
    V3, F2 = cap_delaunay(rings_edge, target_edge)
    # first nphi rows of V3 are the ring itself -> re-use the existing vertices
    remap = np.empty(len(V3), np.int64)
    remap[:nphi] = base_index + np.arange(nphi)
    remap[nphi:] = nv + np.arange(len(V3) - nphi)
    verts.append(V3[nphi:])
    f = remap[F2]
    # orient
    a, b, c = V3[F2[:, 0]], V3[F2[:, 1]], V3[F2[:, 2]]
    n = np.cross(b - a, c - a)
    flip = (n @ outward) < 0
    f[flip] = f[flip][:, ::-1]
    faces.append(f)
    return nv + len(V3) - nphi


def tube_mesh(curve, radii, nphi=96, profile=None, cap_start=True, cap_end=True,
              target_edge=0.8):
    """Watertight (optionally capped) swept tube with well-shaped planar caps."""
    rings = sweep_rings(curve, radii, nphi, profile)
    nsec = len(curve)
    V = rings.reshape(-1, 3)
    verts = [V]
    faces = [loft_faces(nsec, nphi)]
    nv = len(V)
    if cap_start:
        out = curve[0] - curve[1]; out /= np.linalg.norm(out)
        nv = _attach_cap(rings[0], 0, verts, faces, nv, out, target_edge)
    if cap_end:
        out = curve[-1] - curve[-2]; out /= np.linalg.norm(out)
        nv = _attach_cap(rings[-1], (nsec - 1) * nphi, verts, faces, nv, out, target_edge)
    m = trimesh.Trimesh(np.vstack(verts), np.vstack(faces), process=True)
    m.fix_normals()
    return m


def prolate_loft(centre, axis, R, C, base_off, nphi=None, target_edge=0.8):
    """Prolate spheroid truncated by a plane, built as a loft so that the base
    cap and the meridional spacing are both controlled. axis points centre->apex."""
    axis = np.asarray(axis, float); axis /= np.linalg.norm(axis)
    centre = np.asarray(centre, float)
    # sample uniformly along the MERIDIONAL arclength so the apex stays resolved
    t0 = np.arccos(np.clip(-base_off / C, -1, 1))
    tt = np.linspace(t0, 0.0, 6000)
    ds = np.sqrt((C * np.sin(tt)) ** 2 + (R * np.cos(tt)) ** 2)
    sarc = np.r_[0.0, np.cumsum(np.abs(np.diff(tt)) * 0.5 * (ds[:-1] + ds[1:]))]
    n = max(24, int(sarc[-1] / target_edge))
    if nphi is None:
        nphi = max(64, 4 * int(np.ceil(2 * np.pi * R / target_edge / 4)))
    t = np.interp(np.linspace(0, sarc[-1], n), sarc, tt)
    u = C * np.cos(t)
    r = R * np.sin(t)
    r[-1] = 0.0
    curve = centre[None, :] + u[:, None] * axis[None, :]
    rings = sweep_rings(curve, r, nphi)
    verts = [rings.reshape(-1, 3)]
    faces = [loft_faces(n, nphi)]
    nv = n * nphi
    nv = _attach_cap(rings[0], 0, verts, faces, nv, -axis, target_edge)
    m = trimesh.Trimesh(np.vstack(verts), np.vstack(faces), process=True)
    m.merge_vertices()
    m.update_faces(m.nondegenerate_faces())
    m.remove_unreferenced_vertices()
    m.fix_normals()
    return m


def ellipsoid(center, radii, rot=None, subdiv=4):
    """Icosphere-based ellipsoid. rot = 3x3 rotation matrix (columns = local axes)."""
    s = trimesh.creation.icosphere(subdivisions=subdiv, radius=1.0)
    v = s.vertices * np.asarray(radii, float)[None, :]
    if rot is not None:
        v = v @ np.asarray(rot, float).T
    v = v + np.asarray(center, float)[None, :]
    return trimesh.Trimesh(v, s.faces, process=False)


def frame_from_axis(zdir, hint=None):
    """Build an orthonormal 3x3 matrix whose 3rd column is the given axis."""
    z = np.asarray(zdir, float)
    z /= np.linalg.norm(z)
    h = np.array([0.0, 0.0, 1.0]) if hint is None else np.asarray(hint, float)
    if abs(np.dot(h, z)) > 0.95:
        h = np.array([1.0, 0.0, 0.0])
    x = np.cross(h, z)
    x /= np.linalg.norm(x)
    y = np.cross(z, x)
    return np.column_stack([x, y, z])


# ----------------------------------------------------------------------
# booleans
# ----------------------------------------------------------------------
def _clean(m):
    m = m.copy()
    m.merge_vertices()
    m.remove_degenerate_faces() if hasattr(m, 'remove_degenerate_faces') else None
    m.update_faces(m.nondegenerate_faces())
    m.remove_unreferenced_vertices()
    m.fix_normals()
    return m


def union(meshes):
    meshes = [_clean(m) for m in meshes]
    out = trimesh.boolean.union(meshes, engine='manifold')
    return _clean(out)


def difference(a, b):
    out = trimesh.boolean.difference([_clean(a), _clean(b)], engine='manifold')
    return _clean(out)


def intersection(a, b):
    out = trimesh.boolean.intersection([_clean(a), _clean(b)], engine='manifold')
    return _clean(out)


def report(m, name=""):
    print(f"{name:28s} tris={len(m.faces):8d}  watertight={m.is_watertight}  "
          f"winding_ok={m.is_winding_consistent}  vol={m.volume/1000.0:8.2f} mL")
    return m
