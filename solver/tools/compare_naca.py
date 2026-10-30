#!/usr/bin/env python3
"""Compare the NACA66mod cases (C-grid, O-grid, overset) on the y = 0 plane:
surface pressure coefficient along the foil and u / p fields."""
import sys, os, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from readq import read_x, read_q
from scipy.interpolate import griddata
RUN=sys.argv[1] if len(sys.argv)>1 else '..'   # directory holding the case run folders
rho=998.2072; U=1.0; q0=0.5*rho*U*U
import os
cases=[c for c in [('NACA66mod_1phase','C-grid'),('NACA66mod_Ogrid_1phase','O-grid'),('NACA66mod_Overset_1phase','overset')] if os.path.exists(f'{RUN}/{c[0]}/{c[0]}_00000010.q')]
def load(c):
    G=read_x(f'{RUN}/{c}/{c}.x'); Q=read_q(f'{RUN}/{c}/{c}_00000010.q'); return G,Q
data={}
for c,lab in cases:
    G,Q=load(c); pts=[]; wall=None
    for b,(g,q) in enumerate(zip(G,Q)):
        ni,nj,nk=g.shape[1:]
        # y = 0 plane: the first plane of the spanwise direction (k for the foil grids, j for the boxes)
        if g[1][:,0,0].std() < 1e-9 and abs(g[1][:,0,0].mean())<1e-9 and nk==5:   # y varies along k
            x=g[0][:,:,0]; z=g[2][:,:,0]; u=q[0][:,:,0]; p=q[3][:,:,0]; ib=q[-1][:,:,0] if q.shape[0]>4 else np.ones_like(u)
            if b==0: wall=(g[0][:,0,0], g[2][:,0,0], q[3][:,0,0])
        else:                                                   # y along j
            x=g[0][:,0,:]; z=g[2][:,0,:]; u=q[0][:,0,:]; p=q[3][:,0,:]; ib=q[-1][:,0,:] if q.shape[0]>4 else np.ones_like(u)
        m=ib>=0
        pts.append(np.column_stack([x[m],z[m],u[m],p[m]]))
    P=np.vstack(pts)
    # foil position: wall nodes (block 1, j=1; C-grid: i 30..372)
    xw,zw,pw=wall
    if c.startswith('NACA66mod_1phase'): xw,zw,pw=xw[29:372],zw[29:372],pw[29:372]
    zoff=0.5*(zw.max()+zw.min()); xoff=xw.min()
    data[lab]=dict(P=P,xw=xw-xoff,zw=zw-zoff,cp=pw/q0,xoff=xoff,zoff=zoff)
    print(lab,'foil x %.3f..%.3f z %.3f..%.3f -> offsets x %.3f z %.3f, nodes on plane %d'%(xw.min(),xw.max(),zw.min(),zw.max(),xoff,zoff,len(P)))
# ---- surface Cp
fig,ax=plt.subplots(1,2,figsize=(15,5))
for lab,d in data.items():
    ax[0].plot(d['xw'],d['cp'],'.',ms=3,label=lab)
ax[0].set_xlabel('x/c'); ax[0].set_ylabel('Cp'); ax[0].invert_yaxis(); ax[0].legend(); ax[0].set_title('surface pressure coefficient, y = 0, t = 0.02 s')
# ---- profiles on a vertical line x/c = 0.5 and on the wake axis
zs=np.linspace(-0.6,0.6,121); xs=np.linspace(-1,3,161)
for lab,d in data.items():
    P=d['P']; X=P[:,0]-d['xoff']; Z=P[:,1]-d['zoff']
    uz=griddata((X,Z),P[:,2],(0.5*np.ones_like(zs),zs),method='linear')
    ax[1].plot(uz,zs,label=lab+' u(z) at x/c=0.5')
ax[1].set_xlabel('u [m/s]'); ax[1].set_ylabel('z/c'); ax[1].legend(); ax[1].set_title('velocity profile through the foil at x/c = 0.5')
plt.savefig('compare_naca_profiles.png',dpi=80)
# ---- fields
fig,axs=plt.subplots(2,3,figsize=(21,9))
for j,(lab,d) in enumerate(data.items()):
    P=d['P']; X=P[:,0]-d['xoff']; Z=P[:,1]-d['zoff']
    gx,gz=np.meshgrid(np.linspace(-0.5,1.5,201),np.linspace(-0.6,0.6,121))
    for i,(col,ttl,vmin,vmax) in enumerate([(2,'u [m/s]',0.0,1.4),(3,'p [Pa]',-1500,1500)]):
        F=griddata((X,Z),P[:,col],(gx,gz),method='linear')
        if col==3: F=F-np.nanmean(F); ttl='p - <p> [Pa] (level removed)'
        im=axs[i,j].pcolormesh(gx,gz,F,vmin=vmin,vmax=vmax,cmap='jet',shading='auto'); axs[i,j].fill(d['xw'],d['zw'],'w')
        axs[i,j].set_aspect('equal'); axs[i,j].set_title(f'{lab}: {ttl} on y = 0, t = 0.02 s'); plt.colorbar(im,ax=axs[i,j])
plt.savefig('compare_naca_fields.png',dpi=70)
# ---- quantitative: Cp difference on common x/c samples (upper/lower by sign of z)
xs=np.linspace(0.02,0.98,49)
def cp_on(d,side):
    m=(d['zw']>0) if side>0 else (d['zw']<0); o=np.argsort(d['xw'][m]); return np.interp(xs,d['xw'][m][o],d['cp'][m][o])
ref=data['O-grid']
for lab,d in data.items():
    if lab=='O-grid': continue
    for side,nm in [(1,'upper'),(-1,'lower')]:
        dif=cp_on(d,side)-cp_on(ref,side); lvl=dif.mean(); dif0=dif-lvl
        print(f'{lab} vs O-grid, {nm} surface Cp: level offset {lvl:+.3f}; after removing the level: max |dCp| {np.abs(dif0).max():.3f}, rms {np.sqrt((dif0**2).mean()):.3f}  (Cp range of the O-grid {cp_on(ref,side).min():.2f}..{cp_on(ref,side).max():.2f})')
    # velocity profile agreement at x/c = 0.5 and along the wake axis z = 0
    for name,X,Z in [('u(z) at x/c = 0.5',0.5*np.ones_like(zs),zs),('u(x) on z/c = 0.15 (above the foil)',xs,0.15*np.ones_like(xs))]:
        ua=griddata((d['P'][:,0]-d['xoff'],d['P'][:,1]-d['zoff']),d['P'][:,2],(X,Z)); ub=griddata((ref['P'][:,0]-ref['xoff'],ref['P'][:,1]-ref['zoff']),ref['P'][:,2],(X,Z))
        m=np.isfinite(ua)&np.isfinite(ub); print(f'   {name}: max |du| {np.abs(ua-ub)[m].max():.3f} m/s, rms {np.sqrt(((ua-ub)[m]**2).mean()):.3f}')
