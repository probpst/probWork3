#!/usr/bin/env python3
"""morph_wall_frames.py - volume morph of the heart grid for the 97 cycle frames from the delivered
WALL displacements (heart/data/wall_frames/wall_disp_*.npz), with the validated diffusion morpher
(heart/scripts/cycle/diffmorph.py: edge stiffness 1/l^2 (d_ref/d_wall), AMG solve).  Writes
OUT/disp_NNN.npy (global-node displacement, float32, mm) and OUT/frames_log.csv (inverted cells,
min scaled Jacobian of every frame and of the linear midpoint to the previous frame).  Resumable;
--workers N runs N processes on contiguous frame chunks (about 2 GB each).
usage: morph_wall_frames.py --hg HGDIR --cyc CYCDIR --bmesh bmesh.npz --wall WALLDIR --out OUT [--workers 4]
  HGDIR holds assembled_s5.pkl, mb.py, numbering.py (and node_numbering.npz after the first use);
  CYCDIR holds diffmorph.py (its import of numbering must point to HGDIR: run the setup sed first).
"""
import sys, os, csv, time, glob, argparse, pickle, numpy as np

def load_wall(walldir):
    fr, tm, ds = [], [], []
    for f in sorted(glob.glob(os.path.join(walldir, 'wall_disp_*.npz'))):
        z = np.load(f); fr.append(z['frames']); tm.append(z['times']); ds.append(z['disp'])
    fr = np.concatenate(fr); o = np.argsort(fr)
    return fr[o], np.concatenate(tm)[o], np.concatenate(ds)[o]

def work(args, frames):
    sys.path.insert(0, args.hg); sys.path.insert(0, args.cyc)
    from diffmorph import DiffMorph
    A = pickle.load(open(os.path.join(args.hg, 'assembled_s5.pkl'), 'rb')); bm = dict(np.load(args.bmesh))
    fidx, times, W = load_wall(args.wall)
    M = DiffMorph(A, bm['pos'], alpha=1.0)
    if frames[0] == 0:
        np.save(os.path.join(args.out, 'P0.npy'), M.P); np.save(os.path.join(args.out, 'cells.npy'), M.cells)
    log = os.path.join(args.out, f'frames_log_{frames[0]:03d}.csv')
    done = set()
    if os.path.exists(log):
        done = {int(r['frame']) for r in csv.DictReader(open(log))}
    else:
        open(log, 'w').write('frame,time_s,neg,minJ,below005,mid_neg,mid_minJ,sec\n')
    prev = None
    for n in frames:
        f = os.path.join(args.out, f'disp_{n:03d}.npy')
        if n in done and os.path.exists(f):
            prev = np.load(f).astype(np.float64); continue
        t = time.time(); Dn = M.morph(W[n].astype(np.float64)); J = M.jac(M.P + Dn)
        if prev is None and n > 0 and os.path.exists(os.path.join(args.out, f'disp_{n-1:03d}.npy')):
            prev = np.load(os.path.join(args.out, f'disp_{n-1:03d}.npy')).astype(np.float64)
        mneg, mmin = -1, float('nan')
        if prev is not None:
            Jm = M.jac(M.P + 0.5 * (prev + Dn)); mneg = int((Jm <= 0).sum()); mmin = float(Jm.min())
        np.save(f, Dn.astype(np.float32))
        with open(log, 'a') as fo:
            fo.write(f'{n},{times[n]:.6f},{int((J<=0).sum())},{J.min():.5f},{int((J<0.05).sum())},{mneg},{mmin:.5f},{time.time()-t:.0f}\n')
        print(f'[{frames[0]:03d}] frame {n:03d} t={times[n]:.4f} neg {(J<=0).sum()} minJ {J.min():.4f} mid_neg {mneg} ({time.time()-t:.0f}s)', flush=True)
        prev = Dn

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--hg', required=True); ap.add_argument('--cyc', required=True); ap.add_argument('--bmesh', required=True)
    ap.add_argument('--wall', required=True); ap.add_argument('--out', required=True); ap.add_argument('--workers', type=int, default=1)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    fidx, times, W = load_wall(args.wall); nf = len(fidx); del W
    chunks = np.array_split(np.arange(nf), args.workers)
    if args.workers == 1:
        work(args, chunks[0])
    else:
        import multiprocessing as mp
        ps = [mp.Process(target=work, args=(args, c)) for c in chunks]
        for p in ps: p.start()
        for p in ps: p.join()
        if any(p.exitcode != 0 for p in ps): sys.exit('a worker failed')
    # merged log
    rows = []
    for f in sorted(glob.glob(os.path.join(args.out, 'frames_log_*.csv'))):
        rows += list(csv.DictReader(open(f)))
    rows.sort(key=lambda r: int(r['frame']))
    with open(os.path.join(args.out, 'frames_log.csv'), 'w') as fo:
        fo.write('frame,time_s,neg,minJ,below005,mid_neg,mid_minJ,sec\n')
        for r in rows: fo.write(','.join(r[k] for k in ['frame','time_s','neg','minJ','below005','mid_neg','mid_minJ','sec']) + '\n')
    neg = max(int(r['neg']) for r in rows); mneg = max(int(r['mid_neg']) for r in rows)
    print(f'frames {len(rows)} / {nf}: inverted cells max {neg}, min scaled Jacobian {min(float(r["minJ"]) for r in rows):.4f}, '
          f'midpoint inverted max {mneg}')
    if len(rows) == nf and neg == 0 and mneg == 0: print('ALL FRAMES VALID')

if __name__ == '__main__':
    main()
