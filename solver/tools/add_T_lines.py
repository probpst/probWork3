import sys, re
src=open(sys.argv[1]).read().splitlines(); out=[]; nvar=int(sys.argv[3]) if len(sys.argv)>3 else 4
i=0; cnt=0
# a surface record: dir line, 3 range lines, then nvar "flag value" lines; detect by counting after a direction line
while i < len(src):
    l=src[i]; out.append(l)
    toks=l.split('!')[0].split()
    if len(toks)==1 and i+4 < len(src) and len(src[i+1].split('!')[0].split())==2 and len(src[i+2].split('!')[0].split())==2 and len(src[i+3].split('!')[0].split())==2 \
       and abs(int(toks[0]))<=3 and int(toks[0])!=0 and len(src[i+4].split('!')[0].split())>=2:
        # copy ranges + variable lines
        for m in range(1,4+nvar): out.append(src[i+m])
        u=src[i+4].split('!')[0].split()
        if u[0]=='1' and float(u[1])!=0.0: out.append('       1   293.0       ! T: Dirichlet 293 K (inflow)')
        else: out.append('       0   0.0         ! T: zero gradient')
        cnt+=1; i+=4+nvar
    i+=1
open(sys.argv[2],'w').write('\n'.join(out)+'\n'); print('T lines added:',cnt)
