#!/usr/bin/env python3
"""Plot the step history of <Case>_monitor.dat (divergence, iterations, residuals).
Usage: python3 plot_monitor.py Case_monitor.dat"""
import sys, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
d=np.loadtxt(sys.argv[1]); t=d[:,1]
fig,ax=plt.subplots(2,2,figsize=(11,7))
ax[0,0].semilogy(t,d[:,7],label='|div|max'); ax[0,0].semilogy(t,d[:,8],label='|div|avg'); ax[0,0].set_title('divergence [1/s]'); ax[0,0].legend()
ax[0,1].plot(t,d[:,2],label='CFL total'); ax[0,1].plot(t,d[:,3],label='i'); ax[0,1].plot(t,d[:,4],label='j'); ax[0,1].plot(t,d[:,5],label='k'); ax[0,1].set_title('max local CFL'); ax[0,1].legend()
ax[1,0].plot(t,d[:,9],label='momentum sweeps'); ax[1,0].plot(t,d[:,11],label='Poisson iterations'); ax[1,0].set_title('iterations'); ax[1,0].legend()
ax[1,1].semilogy(t,d[:,10],label='res(u*)'); ax[1,1].semilogy(t,d[:,12],label='res(p)'); ax[1,1].semilogy(t,d[:,13],label='upd(p)'); ax[1,1].set_title('residuals'); ax[1,1].legend()
for a in ax.ravel(): a.set_xlabel('t [s]')
plt.tight_layout(); plt.savefig(sys.argv[1].replace('.dat','.png'),dpi=120); print('written',sys.argv[1].replace('.dat','.png'))
