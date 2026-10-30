import sys,time,pickle; sys.path.insert(0,'/home/claude/hg')
import numpy as np, heart, oblayer, mb
dom=sys.argv[1]; M=int(sys.argv[2]); ds=float(sys.argv[3]); gap=float(sys.argv[4]); sc=float(sys.argv[5]) if len(sys.argv)>5 else 1.0; cl=float(sys.argv[6]) if len(sys.argv)>6 else 1.0
H=heart.Build(ds=ds,scale=sc); H.left(); H.right(); H.exclude(gap)
hc={k:oblayer.h_first(*v) for k,v in H.comp_uD.items()}
F=H.FL if dom=='L' else H.FR
core=[b for b in H.blocks if b.dom==dom]; io=[x for x in H.inout if x['dom']==dom]
bad=[b.name for b in core if (mb.hex_jac(b.X)<=0).any()]
print('core invalid:',bad,flush=True)
t=time.time()
bl,I=oblayer.build(core,F,io,M,hc,F.names,clust=cl,n_smooth=40,cosmin=0.3,verbose=True,rfloor=0.05,n_untangle=220)
neg=sum(int((mb.hex_jac(B.X)<=0).sum()) for B in bl)
cells=sum(int(np.prod(np.array(B.dims)-1)) for B in core+bl)
dev=np.abs(I['dev'])
print(f'{dom}: M={M} ds={ds} scale={sc} blocks {len(core)+len(bl)} cells {cells} neg {neg} dev mean {dev.mean():.2f} p95 {np.percentile(dev,95):.2f} max {dev.max():.2f} t={time.time()-t:.0f}s',flush=True)
for B in core+bl:
    for a in ('path','F1','F2','hs','s'): 
        if hasattr(B,a): delattr(B,a)
pickle.dump(dict(blocks=core+bl,inout=io,info={k:I[k] for k in ('L','h','g','which','dev','Wp','pos')},names=F.names,hc=hc,uD=H.comp_uD),open(f'/home/claude/hg/final_{dom}.pkl','wb'))
