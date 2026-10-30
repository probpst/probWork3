#!/usr/bin/env python3
"""Compare the cylinder cases (single grid vs overset) on the z = 0 plane."""
import sys, os, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from readq import read_x, read_q
from scipy.interpolate import griddata
RUN=sys.argv[1] if len(sys.argv)>1 else '..'   # directory holding the case run folders
rho=1000.0; q0=0.5*rho
def plane(g,q,z0=0.0):
    zz=g[2][0,0,:]; k=np.searchsorted(zz,z0)-1; w=(z0-zz[k])/(zz[k+1]-zz[k])
    f=lambda a: (1-w)*a[:,:,k]+w*a[:,:,k+1]
    x=f(g[0]); y=f(g[1]); u=f(q[0]); v=f(q[1]); p=f(q[3]); ib=np.minimum(q[-1][:,:,k],q[-1][:,:,k+1]) if q.shape[0]>5 or (q.shape[0]==5 and 'a' not in 'x') else np.ones_like(u)
    return x,y,u,v,p,ib
def load(c, nvar_base):
    G=read_x(f'{RUN}/{c}/{c}.x'); Q=read_q(f'{RUN}/{c}/{c}_00000010.q')
    pts=[]; cp=[]
    for b,(g,q) in enumerate(zip(G,Q)):
        zz=g[2][0,0,:]; k=np.searchsorted(zz,0.0)-1; w=(0-zz[k])/(zz[k+1]-zz[k]); f=lambda a:(1-w)*a[:,:,k]+w*a[:,:,k+1]
        x=f(g[0]); y=f(g[1]); u=f(q[0]); v=f(q[1]); p=f(q[3])
        ib=np.minimum(q[-1][:,:,k],q[-1][:,:,k+1]) if q.shape[0]>nvar_base else np.ones_like(u)
        m=ib>=0; pts.append(np.column_stack([x[m],y[m],u[m],v[m],p[m]]))
        if b<2: cp.append(np.column_stack([np.degrees(np.arctan2(y[:,0],x[:,0])), p[:,0]/q0]))
    cp=np.vstack(cp); o=np.argsort(cp[:,0])
    return dict(P=np.vstack(pts), th=cp[o,0], cp=cp[o,1])
def run(pair, tag, nvar_base):
    D={lab:load(c,nvar_base) for c,lab in pair if os.path.exists(f'{RUN}/{c}/{c}_00000010.q')}
    if len(D)<2: print(tag,'incomplete'); return
    fig,ax=plt.subplots(1,3,figsize=(21,5))
    for lab,d in D.items(): ax[0].plot(d['th'],d['cp'],'.-',ms=3,label=lab)
    ax[0].set_xlabel('theta [deg] (0 = downstream)'); ax[0].set_ylabel('Cp'); ax[0].legend(); ax[0].set_title(f'{tag}: surface pressure coefficient at z = 0, t = 0.05 s')
    xs=np.linspace(-2.4,6,169); ys=np.linspace(0.55,3.9,68)
    for lab,d in D.items():
        P=d['P']; ax[1].plot(xs,griddata((P[:,0],P[:,1]),P[:,2],(xs,np.zeros_like(xs)),method='linear'),label=lab+' u on y = 0')
        ax[2].plot(griddata((P[:,0],P[:,1]),P[:,2],(np.zeros_like(ys),ys),method='linear'),ys,label=lab+' u on x = 0')
    ax[1].set_xlabel('x [m]'); ax[1].set_ylabel('u [m/s]'); ax[1].legend(); ax[2].set_xlabel('u [m/s]'); ax[2].set_ylabel('y [m]'); ax[2].legend()
    plt.savefig(f'compare_{tag}_profiles.png',dpi=80)
    fig,axs=plt.subplots(2,len(D),figsize=(8*len(D),10))
    gx,gy=np.meshgrid(np.linspace(-2.5,5,151),np.linspace(-3,3,121))
    for j,(lab,d) in enumerate(D.items()):
        P=d['P']
        for i,(col,ttl,vmin,vmax) in enumerate([(2,'u [m/s]',-0.2,1.8),(4,'p [Pa]',-6000,6000)]):
            F=griddata((P[:,0],P[:,1]),P[:,col],(gx,gy),method='linear'); F[np.sqrt(gx**2+gy**2)<0.5]=np.nan
            im=axs[i,j].pcolormesh(gx,gy,F,vmin=vmin,vmax=vmax,cmap='jet',shading='auto'); axs[i,j].set_aspect('equal'); axs[i,j].set_title(f'{lab}: {ttl}, z = 0, t = 0.05 s'); plt.colorbar(im,ax=axs[i,j])
    plt.savefig(f'compare_{tag}_fields.png',dpi=70)
    labs=list(D); a,b=D[labs[0]],D[labs[1]]
    th=np.linspace(-179,179,359); ca=np.interp(th,a['th'],a['cp']); cb=np.interp(th,b['th'],b['cp'])
    print(f'{tag}: {labs[1]} vs {labs[0]}: surface Cp max|dCp| {np.abs(ca-cb).max():.3f} rms {np.sqrt(((ca-cb)**2).mean()):.3f}; Cp range {ca.min():.2f}..{ca.max():.2f}')
    for name,X,Y in [('y=0 axis',xs,np.zeros_like(xs)),('x=0 line',np.zeros_like(ys),ys)]:
        ua=griddata((a['P'][:,0],a['P'][:,1]),a['P'][:,2],(X,Y)); ub=griddata((b['P'][:,0],b['P'][:,1]),b['P'][:,2],(X,Y))
        pa=griddata((a['P'][:,0],a['P'][:,1]),a['P'][:,4],(X,Y)); pb=griddata((b['P'][:,0],b['P'][:,1]),b['P'][:,4],(X,Y))
        m=np.isfinite(ua)&np.isfinite(ub); print(f'   {name}: max|du| {np.abs(ua-ub)[m].max():.3f} m/s, max|dp| {np.abs(pa-pb)[m].max():.0f} Pa (p range {pa[m].min():.0f}..{pa[m].max():.0f})')
run([('Cylinder_1phase','single grid'),('Cylinder_Overset_1phase','overset')],'cylinder_1phase',4)
run([('Cylinder_2phase','single grid'),('Cylinder_Overset_2phase','overset')],'cylinder_2phase',5)
