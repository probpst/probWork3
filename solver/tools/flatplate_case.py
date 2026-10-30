#!/usr/bin/env python3
"""Generate a flat-plate boundary-layer case for probCFD_부산 (update 9 check-up).
   python3 tools/flatplate_case.py <Re_L> <y+_target> <outdir>
   single block: x in [-0.3, 1] (plate 0..1 with bcs flag 11 = wall function),
   y in [0, 0.5] wall normal, z 4 cells spanwise (symmetry planes)."""
import numpy as np, sys, os, shutil
Re=float(sys.argv[1]); yplus=float(sys.argv[2]); out=sys.argv[3]
L=1.0; U=1.0; rho=1000.0; nu=U*L/Re; mu=rho*nu
# first-cell height from the local skin friction at mid plate (Cf = 0.026 Re_x^-1/7)
cf=0.026*(Re*0.5)**(-1/7.); utau=U*np.sqrt(cf/2); y1=yplus*nu/utau*2.0   # cell height = 2 x centre distance
ni,nj,nk=141,49,5
# x: clustered at the leading edge (x=0), cells ahead of the plate coarser
x_up=np.linspace(-0.3,0.0,16)[:-1]
s=np.linspace(0,1,ni-15); beta=2.2; x_pl=np.sinh(beta*s)/np.sinh(beta)
x=np.concatenate([x_up,x_pl])
# y: geometric growth from y1 to H=0.5
H=0.5 if Re>2e5 else 1.5; n=nj-1
def solve_r(y1,H,n):
    lo,hi=1.0001,3.0
    for _ in range(100):
        r=0.5*(lo+hi); tot=y1*(r**n-1)/(r-1)
        if tot>H: hi=r
        else: lo=r
    return 0.5*(lo+hi)
r=solve_r(y1,H,n); y=np.concatenate([[0.0],np.cumsum(y1*r**np.arange(n))]); y[-1]=H
z=np.linspace(0,0.02,nk)
X,Y,Z=np.meshgrid(x,y,z,indexing='ij')
name=os.path.basename(out.rstrip('/')); os.makedirs(out,exist_ok=True)
with open(os.path.join(out,name+'.grd'),'w') as f:
    f.write("%8d\n"%1); f.write("%8d%8d%8d\n"%(ni,nj,nk))
    for a in (X,Y,Z):
        v=a.transpose(2,1,0).ravel(); n4=(len(v)//4)*4
        if n4: np.savetxt(f,v[:n4].reshape(-1,4),fmt='%22.15e')
        if n4<len(v): np.savetxt(f,v[n4:][None,:],fmt='%22.15e')
ip=15+1   # node index where the plate starts (x=0): x_up has 15 points, so node 16 is x=0
bcs=f"""! flat plate, Re_L = {Re:.1e}, first cell y+ ~ {yplus:.0f} (y1 = {y1:.3e} m, growth {r:.3f})
       1                ! number of grid blocks
       7                ! number of boundary surfaces of block 1
       2                ! face direction (j-min): symmetry ahead of the plate
       1  {ip:3d}           ! i-range (nodes)
       1    1           ! j-range (nodes)
       1  {nk:3d}           ! k-range (nodes)
       0   0.0         ! u: zero gradient (slip)
       1   0.000       ! v: Dirichlet 0 (impermeable)
       0   0.0         ! w: zero gradient
       0   0.0         ! p: zero gradient
       2                ! face direction (j-min): PLATE, no-slip wall + wall function
     {ip:3d}  {ni:3d}           ! i-range (nodes)
       1    1           ! j-range (nodes)
       1  {nk:3d}           ! k-range (nodes)
      11   0.000       ! u: flag 11 = no-slip wall with wall function
      11   0.000       ! v:
      11   0.000       ! w:
       0   0.0         ! p: zero gradient
      -2                ! face direction (j-max): top, slip
       1  {ni:3d}           ! i-range (nodes)
     {nj:3d}  {nj:3d}           ! j-range (nodes)
       1  {nk:3d}           ! k-range (nodes)
       0   0.0         ! u: zero gradient
       1   0.000       ! v: Dirichlet 0
       0   0.0         ! w: zero gradient
       0   0.0         ! p: zero gradient
       1                ! face direction (i-min): inlet
       1    1           ! i-range (nodes)
       1  {nj:3d}           ! j-range (nodes)
       1  {nk:3d}           ! k-range (nodes)
       1   1.000       ! u: Dirichlet 1 m/s
       1   0.000       ! v: Dirichlet 0
       1   0.000       ! w: Dirichlet 0
       0   0.0         ! p: zero gradient
      -1                ! face direction (i-max): outlet
     {ni:3d}  {ni:3d}           ! i-range (nodes)
       1  {nj:3d}           ! j-range (nodes)
       1  {nk:3d}           ! k-range (nodes)
       2   0.0         ! u: outflow (zero gradient, mass-flux corrected)
       2   0.0         ! v: outflow
       2   0.0         ! w: outflow
       1   0.000       ! p: Dirichlet 0
       3                ! face direction (k-min): symmetry
       1  {ni:3d}           ! i-range (nodes)
       1  {nj:3d}           ! j-range (nodes)
       1    1           ! k-range (nodes)
       0   0.0         ! u: zero gradient
       0   0.0         ! v: zero gradient
       1   0.000       ! w: Dirichlet 0
       0   0.0         ! p: zero gradient
      -3                ! face direction (k-max): symmetry
       1  {ni:3d}           ! i-range (nodes)
       1  {nj:3d}           ! j-range (nodes)
     {nk:3d}  {nk:3d}           ! k-range (nodes)
       0   0.0         ! u: zero gradient
       0   0.0         ! v: zero gradient
       1   0.000       ! w: Dirichlet 0
       0   0.0         ! p: zero gradient
"""
open(os.path.join(out,name+'.bcs'),'w').write(bcs)
here=os.path.dirname(os.path.abspath(__file__)); base=os.path.join(here,'..')
nml=open(os.path.join(base,'NACA66mod_Ogrid_1phase','probCFD.nml')).read()
import re
def setv(s,key,val):
    return re.sub(r'^(\s*%s\s*=\s*)[^!\n]*'%key, lambda m: m.group(1)+str(val).ljust(12), s, flags=re.M)
nml=setv(nml,'Case_name','"%s"'%name); nml=setv(nml,'Del_T','0.01'); nml=setv(nml,'End_T','6.0')
nml=setv(nml,'mu_w','%.6e'%mu); nml=setv(nml,'Output_Intvl','100'); nml=setv(nml,'TURB_model','4')
nml=setv(nml,'t_ramp','0.05'); nml=setv(nml,'Cnvt_mode','2'); nml=setv(nml,'Log_Intvl','10')
open(os.path.join(out,'probCFD.nml'),'w').write(nml)
shutil.copy(os.path.join(base,'NACA66mod_Ogrid_1phase','Input_polyMesh_probCFD.csv'),os.path.join(out,'Input_polyMesh_probCFD.csv'))
ittc=0.075/(np.log10(Re)-2)**2
open(os.path.join(out,'README.md'),'w').write(f"""# {name}: flat-plate boundary layer, Re_L = {Re:.1e}

Single block {ni} x {nj} x {nk} nodes, plate on y = 0 for 0 <= x <= 1 m (bcs flag 11:
no-slip wall + wall function), symmetry ahead of the plate, slip top at y = {H} m,
uniform inlet U = 1 m/s at x = -0.3 m, outlet at x = 1 m, spanwise symmetry planes.
Water rho = 1000 kg/m3, mu = {mu:.3e} Pa s (nu = {nu:.3e} m2/s).  First cell height
{y1:.3e} m (centre at y+ ~ {yplus:.0f} for the mid-plate skin friction), growth ratio {r:.3f}.
Model: k-omega SST with the wall function (TURB_model = 4).
Reference: ITTC-1957 CF = 0.075 / (log10 Re - 2)^2 = {ittc:.5f}.
The monitor line "wall-function y+ min/max ... mean Cf" gives the plate-averaged CF.
""")
print(name,'Re %.1e nu %.3e y1 %.3e r %.3f ITTC CF %.5f'%(Re,nu,y1,r,ittc))
