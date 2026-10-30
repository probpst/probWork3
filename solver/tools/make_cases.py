#!/usr/bin/env python3
"""Builds the 7 tutorial case folders (bcs, csv) of probCFD_Busan update 4 from the
supplied grids.  Usage: make_cases.py <dir with the unzipped Cylinder/, Cylinder_Overset/,
NACA66mod/, NACA66mod_Ogrid/, NACA66mod_Overset/ folders> <output dir>"""
import os, shutil, sys
T=sys.argv[1]; OUT=sys.argv[2]
def rec(d, ir, jr, kr, lines, title):
    n={1:'i-min',-1:'i-max',2:'j-min',-2:'j-max',3:'k-min',-3:'k-max'}[d]
    L=[f"  ! --- {n}: {title}", f"      {d:2d}                ! face direction ({n})",
       f"    {ir[0]:4d} {ir[1]:4d}           ! i-range (nodes)",
       f"    {jr[0]:4d} {jr[1]:4d}           ! j-range (nodes)",
       f"    {kr[0]:4d} {kr[1]:4d}           ! k-range (nodes)"]
    for fl,val,c in lines: L.append(f"       {fl}   {val:<12}! {c}")
    return L
# ---- variable line sets ---------------------------------------------------
def wall(tp):   return [('1','0.000','u: Dirichlet 0 (no-slip wall)'),('1','0.000','v: Dirichlet 0'),('1','0.000','w: Dirichlet 0'),('0','0.0','p: zero gradient')]+([('0','0.0','alpha: zero gradient')] if tp else [])
def inlet(tp):  return [('1','1.000','u: Dirichlet 1 m/s (free stream along +x)'),('1','0.000','v: Dirichlet 0'),('1','0.000','w: Dirichlet 0'),('0','0.0','p: zero gradient')]+([('2','0.0','alpha: still-water profile (liquid below FreeSurf_Z)')] if tp else [])
def outlet1():  return [('0','0.0','u: zero gradient'),('0','0.0','v: zero gradient'),('0','0.0','w: zero gradient'),('1','0.000','p: Dirichlet 0 Pa (anchors the pressure level)')]
def outlet2():  return [('2','0.0','u: outflow (zero gradient + global mass-flux correction)'),('2','0.0','v: outflow'),('2','0.0','w: outflow'),('0','0.0','p: zero gradient'),('0','0.0','alpha: zero gradient')]
def slip(ax,tp,c='free slip'):
    L=[('0','0.0',f'u: zero gradient ({c})'),('0','0.0',f'v: zero gradient ({c})'),('0','0.0',f'w: zero gradient ({c})'),('0','0.0','p: zero gradient')]
    L[ax]=('1','0.000',['u','v','w'][ax]+': Dirichlet 0 (impermeable)')
    return L+([('0','0.0','alpha: zero gradient')] if tp else [])
def topopen():  return [('0','0.0','u: zero gradient (open top, air)'),('0','0.0','v: zero gradient'),('0','0.0','w: zero gradient'),('1','0.000','p: Dirichlet 0 Pa (reference pressure)'),('0','0.0','alpha: zero gradient')]
def sym(ax):    return slip(ax,False,'symmetry')
def write(case, header, blocks, gridsrc, csv):
    d=f'{OUT}/{case}'; os.makedirs(d, exist_ok=True)
    shutil.copy(gridsrc, f'{d}/{case}.grd')
    # strip CR
    g=open(f'{d}/{case}.grd','rb').read().replace(b'\r',b''); open(f'{d}/{case}.grd','wb').write(g)
    L=[f"! probCFD boundary conditions for case '{case}' (grid {case}.grd)"]+header+[f"       {len(blocks)}                ! number of grid blocks"]
    for b,(btitle,recs) in enumerate(blocks):
        L.append(f"! ==================== block {b+1}: {btitle}")
        L.append(f"       {len(recs)}                ! number of boundary surfaces of block {b+1}")
        for r in recs: L+=r
    open(f'{d}/{case}.bcs','w').write("\n".join(L)+"\n")
    open(f'{d}/Input_polyMesh_probCFD.csv','w').write(csv)
    print('wrote', case, sum(len(r) for _,r in blocks), 'surfaces')

LAY="""! Record layout per face: direction (|d| = 1 i, 2 j, 3 k; + min, - max), i-, j-,
!   k-range (NODE indices), then one 'flag value' line per variable u, v, w, p
!   [, alpha].  flag 1 = Dirichlet (velocity in global Cartesian m/s, p in Pa),
!   flag 0 = zero gradient, velocity flag 2 = outflow (zero gradient + global
!   mass-flux correction), alpha flag 2 = still-water profile at an inlet.
! Block-to-block interfaces are NOT listed (the solver detects them, .bcs_b2b);
! faces neither listed nor interfaces are OVERSET boundaries (overset cases)."""
# ======================= Cylinder (4-block single grid) =====================
CYL_H=["! Surface-piercing circular cylinder D = 1 m (axis = z): O-grid halves (blocks 1, 2,",
"! j = 1 wall, j = 55 at r = 1.5 m) + two outer blocks 3 (y > 0) and 4 (y < 0);",
"! domain x -2.5 .. 9.5, y -4 .. 4, z -2 .. 1.5 m; k is the vertical direction.",
"! Interfaces: blocks 1/2 i-min, i-max (y = 0 planes), j-max (r = 1.5); blocks 3/4 j-min.", LAY]
def cyl_blocks(tp):
    top = topopen() if tp else slip(2,tp,'free-slip lid')
    b12=lambda: [rec(2,(1,20),(1,1),(1,35),wall(tp),'cylinder surface, no-slip wall'),
                 rec(3,(1,20),(1,55),(1,1),slip(2,tp),'bottom z = -2 m'),
                 rec(-3,(1,20),(1,55),(35,35),top,'top z = 1.5 m'+(' (air, open)' if tp else ' (single phase: free-slip lid)'))]
    b3=[rec(1,(1,1),(1,30),(1,35),inlet(tp),'inlet x = -2.5 m'),
        rec(-1,(62,62),(1,30),(1,35),outlet2() if tp else outlet1(),'outlet x = 9.5 m'),
        rec(-2,(1,62),(30,30),(1,35),slip(1,tp),'side y = +4 m'),
        rec(3,(1,62),(1,30),(1,1),slip(2,tp),'bottom'), rec(-3,(1,62),(1,30),(35,35),top,'top')]
    b4=[rec(1,(1,1),(1,30),(1,35),outlet2() if tp else outlet1(),'outlet x = 9.5 m (i reversed in block 4)'),
        rec(-1,(62,62),(1,30),(1,35),inlet(tp),'inlet x = -2.5 m'),
        rec(-2,(1,62),(30,30),(1,35),slip(1,tp),'side y = -4 m'),
        rec(3,(1,62),(1,30),(1,1),slip(2,tp),'bottom'), rec(-3,(1,62),(1,30),(35,35),top,'top')]
    return [('O-grid half y >= 0, 20 x 55 x 35 nodes',b12()),('O-grid half y <= 0, 20 x 55 x 35 nodes',b12()),
            ('outer block y >= 0, 62 x 30 x 35 nodes',b3),('outer block y <= 0, 62 x 30 x 35 nodes (i runs from x = 9.5 to -2.5)',b4)]
CSV4="probCFD polyMesh Input,block-1,block-2,block-3,block-4\ni_split,1,1,1,1\nj_split,1,1,1,1\nk_split,2,2,2,2\nn_of_BLK,2,2,2,2\n"
write('Cylinder_1phase', CYL_H+['! SINGLE PHASE (water).'], cyl_blocks(False), T+'/Cylinder/Cylinder.grd', CSV4)
write('Cylinder_2phase', CYL_H+['! TWO PHASE: water below z = 0 (FreeSurf_Z of the name list), air above.'], cyl_blocks(True), T+'/Cylinder/Cylinder.grd', CSV4)
# ======================= Cylinder overset (3 blocks) =========================
OV_H=["! Surface-piercing circular cylinder D = 1 m (axis = z), OVERSET: O-grid halves",
"! (blocks 1, 2: j = 1 wall, outer boundary r = 1.57 m NOT listed = overset boundary)",
"! inside the Cartesian background block 3 (x -2.5 .. 9.5, y -4 .. 4, z -2 .. 1.5 m,",
"! cells 0.13 x 0.14 x 0.10 m); k is the vertical direction.",
"! Interfaces: blocks 1/2 i-min, i-max (y = 0 planes).", LAY]
def ov_blocks(tp):
    top = topopen() if tp else slip(2,tp,'free-slip lid')
    b12=lambda: [rec(2,(1,20),(1,1),(1,35),wall(tp),'cylinder surface, no-slip wall'),
                 rec(3,(1,20),(1,55),(1,1),slip(2,tp),'bottom z = -2 m'),
                 rec(-3,(1,20),(1,55),(35,35),top,'top z = 1.5 m'+(' (air, open)' if tp else ' (single phase: free-slip lid)'))]
    b3=[rec(1,(1,1),(1,59),(1,35),inlet(tp),'inlet x = -2.5 m'),
        rec(-1,(62,62),(1,59),(1,35),outlet2() if tp else outlet1(),'outlet x = 9.5 m'),
        rec(2,(1,62),(1,1),(1,35),slip(1,tp),'side y = -4 m'), rec(-2,(1,62),(59,59),(1,35),slip(1,tp),'side y = +4 m'),
        rec(3,(1,62),(1,59),(1,1),slip(2,tp),'bottom'), rec(-3,(1,62),(1,59),(35,35),top,'top')]
    return [('O-grid half y >= 0, 20 x 55 x 35 nodes',b12()),('O-grid half y <= 0, 20 x 55 x 35 nodes',b12()),
            ('background, 62 x 59 x 35 nodes',b3)]
CSV3="probCFD polyMesh Input,block-1,block-2,block-3\ni_split,1,1,2\nj_split,1,1,1\nk_split,2,2,2\nn_of_BLK,2,2,4\n"
write('Cylinder_Overset_1phase', OV_H+['! SINGLE PHASE (water).'], ov_blocks(False), T+'/Cylinder_Overset/Cylinder_Overset.grd', CSV3)
write('Cylinder_Overset_2phase', OV_H+['! TWO PHASE: water below z = 0 (FreeSurf_Z), air above.'], ov_blocks(True), T+'/Cylinder_Overset/Cylinder_Overset.grd', CSV3)
# ======================= NACA66mod C-grid (2 blocks) ==========================
C_H=["! NACA66(mod) section, chord 1 m along x, span y 0 .. 0.2 m (5 nodes), SINGLE PHASE.",
"! Block 1: C-grid 401 x 85 x 5; i from the outlet plane x = 21 (upper wake cut, i = 1)",
"!   around the section to x = 21 (lower wake cut, i = 401); j = 1 = foil surface for",
"!   i = 30 .. 372, wake cut lines for i = 1..30 / 372..401 (interfaces with block 2);",
"!   j = 85 = outer C boundary (~20 chords); k spanwise.",
"! Block 2: wake slit 4 x 30 x 5 behind the blunt trailing edge (i across the slit,",
"!   j along x from 1 to 21, k spanwise); j = 1 = trailing edge, i = 1 / i = 4 = interfaces.", LAY]
C_H+=["! OUTLETS (x = 21): pressure outlets p = 0 with zero-gradient velocity (like the O-grid and",
"!   overset cases); the far-field top/bottom segments carry the free stream."]
b1=[rec(1,(1,1),(1,85),(1,5),outlet1(),'outlet x = 21 (upper half)'), rec(-1,(401,401),(1,85),(1,5),outlet1(),'outlet x = 21 (lower half)'),
    rec(2,(30,372),(1,1),(1,5),wall(False),'foil surface (i = 30 .. 372), no-slip wall'),
    rec(-2,(29,374),(85,85),(1,5),inlet(False),'outer C boundary (i = 29 .. 374): inlet'),
    rec(-2,(1,29),(85,85),(1,5),inlet(False),'top far-field segment z = 21.8 (x 1.2 .. 21): free stream (1, 0, 0)'),
    rec(-2,(374,401),(85,85),(1,5),inlet(False),'bottom far-field segment z = -18.8 (x 1.4 .. 21): free stream (1, 0, 0)'),
    rec(3,(1,401),(1,85),(1,1),sym(1),'y = 0 symmetry'), rec(-3,(1,401),(1,85),(5,5),sym(1),'y = 0.2 symmetry')]
b2=[rec(2,(1,4),(1,1),(1,5),wall(False),'blunt trailing edge, no-slip wall'), rec(-2,(1,4),(30,30),(1,5),outlet1(),'outlet x = 21'),
    rec(3,(1,4),(1,30),(1,1),sym(1),'y = 0 symmetry'), rec(-3,(1,4),(1,30),(5,5),sym(1),'y = 0.2 symmetry')]
write('NACA66mod_1phase', C_H, [('C-grid, 401 x 85 x 5 nodes',b1),('wake slit, 4 x 30 x 5 nodes',b2)], T+'/NACA66mod/NACA66mod.grd',
      "probCFD polyMesh Input,block-1,block-2\ni_split,2,1\nj_split,1,1\nk_split,1,1\nn_of_BLK,2,1\n")
# ======================= NACA66mod O-grid (1 block) ===========================
O_H=["! NACA66(mod) section, chord 1 m along x, span y 0 .. 0.2 m (5 nodes), SINGLE PHASE.",
"! O-grid 344 x 81 x 5: i around the section (i = 1 / 344 = O-cut interface, auto-",
"! detected), j = 1 foil surface, j = 81 outer boundary at ~14 chords: the upstream",
"! half (i = 170 .. 344 and the first cell) is the inlet, the downstream half the outlet.", LAY]
bo=[rec(2,(1,344),(1,1),(1,5),wall(False),'foil surface, no-slip wall'),
    rec(-2,(2,170),(81,81),(1,5),outlet1(),'outer boundary, downstream half (cells 2 .. 169): outlet'),
    rec(-2,(170,344),(81,81),(1,5),inlet(False),'outer boundary, upstream half (cells 170 .. 343): inlet'),
    rec(-2,(1,2),(81,81),(1,5),inlet(False),'outer boundary, cell 1 (x = 0.2, z = -13, downstream of the cut): inlet'),
    rec(3,(1,344),(1,81),(1,1),sym(1),'y = 0 symmetry'), rec(-3,(1,344),(1,81),(5,5),sym(1),'y = 0.2 symmetry')]
write('NACA66mod_Ogrid_1phase', O_H, [('O-grid, 344 x 81 x 5 nodes',bo)], T+'/NACA66mod_Ogrid/NACA66mod_Ogrid.grd',
      "probCFD polyMesh Input,block-1\ni_split,2\nj_split,1\nk_split,1\nn_of_BLK,2\n")
# ======================= NACA66mod overset (3 blocks) =========================
V_H=["! NACA66(mod) section, chord 1 m along x, span y 0 .. 0.2 m (5 nodes), SINGLE PHASE,",
"! OVERSET: block 1 foil O-grid 344 x 44 x 5 (j = 1 wall, i = 1 / 344 O-cut interface,",
"! outer boundary j = 44 NOT listed = overset), block 2 refinement box 120 x 5 x 120",
"! (i = x, j = y, k = z; only its y planes are listed), block 3 background 241 x 5 x 80",
"! (x -3 .. 10.5, z -2.73 .. 2.05 m).", LAY]
v1=[rec(2,(1,344),(1,1),(1,5),wall(False),'foil surface, no-slip wall'), rec(3,(1,344),(1,44),(1,1),sym(1),'y = 0 symmetry'), rec(-3,(1,344),(1,44),(5,5),sym(1),'y = 0.2 symmetry')]
v2=[rec(2,(1,120),(1,1),(1,120),sym(1),'y = 0 symmetry'), rec(-2,(1,120),(5,5),(1,120),sym(1),'y = 0.2 symmetry')]
v3=[rec(1,(1,1),(1,5),(1,80),inlet(False),'inlet x = -3'), rec(-1,(241,241),(1,5),(1,80),outlet1(),'outlet x = 10.5'),
    rec(2,(1,241),(1,1),(1,80),sym(1),'y = 0 symmetry'), rec(-2,(1,241),(5,5),(1,80),sym(1),'y = 0.2 symmetry'),
    rec(3,(1,241),(1,5),(1,1),slip(2,False),'bottom z = -2.73'), rec(-3,(1,241),(1,5),(80,80),slip(2,False),'top z = 2.05')]
write('NACA66mod_Overset_1phase', V_H, [('foil O-grid, 344 x 44 x 5 nodes',v1),('refinement box, 120 x 5 x 120 nodes',v2),('background, 241 x 5 x 80 nodes',v3)],
      T+'/NACA66mod_Overset/NACA66mod_Overset.grd', "probCFD polyMesh Input,block-1,block-2,block-3\ni_split,2,2,2\nj_split,1,1,1\nk_split,1,1,1\nn_of_BLK,2,2,2\n")
