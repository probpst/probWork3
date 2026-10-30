#!/usr/bin/env python3
"""u, w, p contours on the mid-span symmetry plane (y = 0.1 C-grid/overset, y = 0.5 O-grid)
of the NACA66mod cases.  usage: plot_mid.py <step> <dir1:label1> <dir2:label2> ... (dirs hold <case>.x/.grd + .q)"""
import sys, os, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0,'/home/claude/probCFD_ovst/tools'); from readq import read_x, read_q
from scipy.interpolate import griddata
def read_grd(fn):
    t=open(fn).read().split(); nb=int(t[0]); dims=np.array(t[1:1+3*nb],int).reshape(nb,3); p=1+3*nb; out=[]
    for b in range(nb):
        ni,nj,nk=dims[b]; n=ni*nj*nk; a=np.array(t[p:p+3*n],float).reshape(3,nk,nj,ni).transpose(0,3,2,1); p+=3*n; out.append(a)
    return out
step=int(sys.argv[1]); items=[a.split(':') for a in sys.argv[2:]]
rho=998.2072; q0=0.5*rho; out=sys.argv[0].replace('plot_mid.py','')
fig,axs=plt.subplots(3,len(items),figsize=(7*len(items),15)); fig2,ax2=plt.subplots(1,2,figsize=(16,6))
gx,gz=np.meshgrid(np.linspace(-0.5,1.5,301),np.linspace(-0.6,0.6,181))
for j,(d,lab) in enumerate(items):
    case=[f[:-2] for f in os.listdir(d) if f.endswith('.x')]
    case=case[0] if case else [f[:-4] for f in os.listdir(d) if f.endswith('.grd')][0]
    G=read_x(f'{d}/{case}.x') if os.path.exists(f'{d}/{case}.x') else read_grd(f'{d}/{case}.grd')
    Q=read_q(f'{d}/{case}_{step:08d}.q'); pts=[]
    for b,(g,q) in enumerate(zip(G,Q)):
        nk=g.shape[3]; m=2                                  # mid plane index (3rd of 5 planes)
        if nk==5: x=g[0][:,:,m]; z=g[2][:,:,m]; u=q[0][:,:,m]; w=q[2][:,:,m]; p=q[3][:,:,m]; ib=q[-1][:,:,m] if q.shape[0]>4 else np.ones_like(u); y=g[1][0,0,m]
        else:     x=g[0][:,m,:]; z=g[2][:,m,:]; u=q[0][:,m,:]; w=q[2][:,m,:]; p=q[3][:,m,:]; ib=q[-1][:,m,:] if q.shape[0]>4 else np.ones_like(u); y=g[1][0,m,0]
        if b==0: xw=g[0][:,0,m]; zw=g[2][:,0,m]; pw=q[3][:,0,m]
        msk=ib>=0; pts.append(np.column_stack([x[msk],z[msk],u[msk],w[msk],p[msk]]))
    P=np.vstack(pts)
    if 'C-grid' in lab or case.startswith('NACA66mod_1phase'): xw,zw,pw=xw[29:372],zw[29:372],pw[29:372]
    zoff=0.5*(zw.max()+zw.min()); xoff=xw.min(); X=P[:,0]-xoff; Z=P[:,1]-zoff
    lvl=np.nanmean(griddata((X,Z),P[:,4],(gx,gz),method='linear'))
    for i,(col,ttl,vmin,vmax) in enumerate([(2,'u [m/s]',0.0,1.4),(3,'w [m/s]',-0.4,0.4),(4,'p [Pa]',-600,600)]):
        F=griddata((X,Z),P[:,col],(gx,gz),method='linear')
        im=axs[i,j].contourf(gx,gz,F,levels=np.linspace(vmin,vmax,25),cmap='jet',extend='both'); axs[i,j].fill(xw-xoff,zw-zoff,'w',zorder=3)
        axs[i,j].set_aspect('equal'); axs[i,j].set_title(f'{lab}: {ttl} on y = {y:.2f}, t = {step*0.002:.3f} s'); plt.colorbar(im,ax=axs[i,j])
    ax2[0].plot(xw-xoff,pw/q0,'.',ms=2,label=lab)
    zs=np.linspace(-0.6,0.6,241); ax2[1].plot(griddata((X,Z),P[:,2],(0.5*np.ones_like(zs),zs)),zs,label=lab)
    print(f'{lab:26s} y = {y:.2f}: u {P[:,2].min():.3f}..{P[:,2].max():.3f}  w {P[:,3].min():.3f}..{P[:,3].max():.3f}  p {P[:,4].min():.0f}..{P[:,4].max():.0f} Pa (mean over the window {lvl:.0f})  wall Cp {pw.min()/q0:.2f}..{pw.max()/q0:.2f}')
ax2[0].invert_yaxis(); ax2[0].legend(); ax2[0].set_xlabel('x/c'); ax2[0].set_ylabel('Cp'); ax2[0].set_title(f'surface pressure coefficient (absolute), t = {step*0.002:.3f} s'); ax2[0].grid(alpha=.3)
ax2[1].legend(); ax2[1].set_xlabel('u [m/s]'); ax2[1].set_ylabel('z/c'); ax2[1].set_title('u(z) at x/c = 0.5'); ax2[1].grid(alpha=.3)
fig.savefig(f'{out}naca_midplane_fields_step{step}.png',dpi=65); fig2.savefig(f'{out}naca_midplane_profiles_step{step}.png',dpi=75)
