"""Greedy frame selection: linear interpolation between consecutive frames must reproduce every
sampled wall position within beta*max(h_first, h_floor); frame spacing <= dt_max."""
import numpy as np,sys,json
D=np.load('Dsamples.npy',mmap_mode='r'); NS=D.shape[0]; h=np.load('bmesh.npz')['h']
T_cycle=60/70.0
def pick(beta,hf,dtmax_ms):
    tol=beta*np.maximum(h,hf); mmax=max(1,int(dtmax_ms/1000/T_cycle*NS)); fr=[0]; i=0
    while i<NS:
        best=i+1; j=i+2
        while j<=min(i+mmax,NS):
            Dj=D[j%NS]; Di=D[i]; ok=True
            for k in range(i+1,j):
                w=(k-i)/(j-i); err=np.linalg.norm(D[k]-((1-w)*Di+w*Dj),axis=1)
                if (err>tol).any(): ok=False; break
            if not ok: break
            best=j; j+=1
        i=best
        if i<NS: fr.append(i)
    return fr
res={}
for beta,hf,dt in [(0.2,0.15,10),(0.1,0.15,10),(0.2,0.25,10)]:
    fr=pick(beta,hf,dt); res[f'{beta}_{hf}_{dt}']=fr
    gaps=np.diff(fr+[NS])*T_cycle/NS*1000
    print(f'beta={beta} h_floor={hf} dt_max={dt} ms -> {len(fr)} frames, spacing {gaps.min():.1f}-{gaps.max():.1f} ms',flush=True)
json.dump(res,open('frame_sets.json','w'))
