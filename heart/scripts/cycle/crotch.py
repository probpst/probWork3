import numpy as np,pickle,sys; sys.path.insert(0,'/home/claude/hg'); import heart
def crotch_points(A):
    B={b.name:b for b in A['blocks']}; H=heart.RES; sc=0.85
    Nm1=max(4,int(round(H['Nm1']*sc))); Nrv=max(4,int(round(H['Nrvot']*sc)))
    G=H.get('Ngap',0); pts=[B['LV'].X[Nm1+G//2,:,0], B['RV'].X[Nrv+H.get('Ngap_rv',0)//2,:,0]]          # split lines on the LV / RV base faces (kmin)
    PA=B['PA'].X[-1]; pts+= [PA[PA.shape[0]//2,:],PA[:,PA.shape[1]//2]]   # PA end-face mid lines (RPA|LPA split)
    doms=['L','R','R','R']
    import os
    if os.environ.get('ISTHMUS','1')=='1':
        TV=B['TV'].X; pts+=[TV[:14,:,0].reshape(-1,3)]; doms+=['R']       # tricuspid lower side near the RA (isthmus ridge)
    P=np.vstack(pts); D=np.concatenate([[d]*len(p) for p,d in zip(pts,doms)])
    return (P,D)

def crotch_groups(A):
    B={b.name:b for b in A['blocks']}; H=heart.RES; sc=0.85
    Nm1=max(4,int(round(H['Nm1']*sc))); Nrv=max(4,int(round(H['Nrvot']*sc)))
    G=H.get('Ngap',0); return [B['LV'].X[Nm1+G//2,:,0], B['RV'].X[Nrv+H.get('Ngap_rv',0)//2,:,0]]
