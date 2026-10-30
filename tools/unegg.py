import struct,sys,os,zlib,bz2,lzma,glob,datetime
vols=sorted(glob.glob(os.path.join(sys.argv[1],'Heart.vol*.egg')),key=lambda p:int(p.split('.vol')[1].split('.')[0]))
out=sys.argv[2]
stream=bytearray()
for v in vols:
    b=open(v,'rb').read()
    assert b[:4]==b'EGGA'
    p=14
    while True:
        magic=struct.unpack_from('<I',b,p)[0]
        if magic==0x08E28222: p+=4; break
        flag=b[p+4]; 
        if flag&1: size=struct.unpack_from('<I',b,p+5)[0]; p+=9+size
        else: size=struct.unpack_from('<H',b,p+5)[0]; p+=7+size
    stream+=b[p:]
print('stream bytes',len(stream))
b=bytes(stream); p=0
def extras(p):
    ex={}
    while True:
        magic=struct.unpack_from('<I',b,p)[0]
        if magic==0x08E28222: return ex,p+4
        flag=b[p+4]
        if flag&1: size=struct.unpack_from('<I',b,p+5)[0]; d=b[p+9:p+9+size]; p+=9+size
        else: size=struct.unpack_from('<H',b,p+5)[0]; d=b[p+7:p+7+size]; p+=7+size
        ex[magic]=(flag,d)
files=[]
while p<len(b):
    magic=struct.unpack_from('<I',b,p)[0]
    if magic==0x0A8590E3:
        fid,flen=struct.unpack_from('<IQ',b,p+4); ex,p=extras(p+16)
        name=None; mtime=None; attr=None
        if 0x0A8591AC in ex:
            flag,d=ex[0x0A8591AC]; q=0
            if flag&8: q+=2
            if flag&4: q+=4
            name=d[q:].decode('utf-8',errors='replace')
        if 0x2C86950B in ex:
            ft,attr=struct.unpack_from('<QB',ex[0x2C86950B][1]); mtime=datetime.datetime(1601,1,1)+datetime.timedelta(microseconds=ft/10)
        cur={'name':name,'len':flen,'attr':attr,'mtime':mtime,'blocks':[],'enc':0x08D1470F in ex}
        files.append(cur)
    elif magic==0x02B50C13:
        meth,hint,usz,csz,crc=struct.unpack_from('<BBIII',b,p+4); ex,p=extras(p+18)
        data=b[p:p+csz]; p+=csz
        cur['blocks'].append((meth,usz,csz,crc,data))
    elif magic==0x08E28222:
        p+=4; break
    else:
        print('unknown magic %08x at %d'%(magic,p)); break
print('files',len(files))
for f in files:
    path=os.path.join(out,f['name'])
    if f['attr'] is not None and f['attr']&0x10 or f['len']==0 and not f['blocks'] and f['name'].endswith('/'):
        os.makedirs(path,exist_ok=True); print('dir ',f['name']); continue
    os.makedirs(os.path.dirname(path) or out,exist_ok=True)
    with open(path,'wb') as fo:
        for meth,usz,csz,crc,data in f['blocks']:
            if meth==0: d=data
            elif meth==1: d=zlib.decompress(data,-15)
            elif meth==2: d=bz2.decompress(data)
            elif meth==4: d=lzma.decompress(data,format=lzma.FORMAT_RAW,filters=[{'id':lzma.FILTER_LZMA1}]) if False else lzma.LZMADecompressor(format=lzma.FORMAT_ALONE).decompress(data)
            else: raise SystemExit('unsupported method %d for %s'%(meth,f['name']))
            assert len(d)==usz and (zlib.crc32(d)&0xffffffff)==crc, ('bad block',f['name'],len(d),usz)
            fo.write(d)
    print('%10d  %s  %s  enc=%s meth=%s'%(f['len'],f['mtime'],f['name'],f['enc'],sorted({bl[0] for bl in f['blocks']})))
