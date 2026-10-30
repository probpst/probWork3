#!/usr/bin/env python3
"""write_frames_bin.py - write a probCFD frame file (morph_frames.F90, 'FRAM' format) from
per-frame node displacements.

Two input forms:
  (a) heart: --global DIR  : DIR/disp_NNN.npy = displacement of the GLOBAL nodes (float32, N x 3, mm)
      with the block<->global map node_numbering.npz (gid, inv) and the block dims from heart_grid.g;
      frame times from cycle_frames.csv (time_s column); periodic with period 60/bpm.
  (b) generic: --blocks DIR : DIR/frame_NNN.npz holding one array per block 'b001'..(3, ni, nj, nk)
      and 'times.npy' (seconds) in DIR.
File layout (little-endian): int32 magic 0x4652414D, int32 version 1, int32 nblk, int32 nframes,
int32 periodic, float64 period [s], float64 scale [m per file unit], int32 dims(3, nblk) (i, j, k per
block, block-major), float64 times(nframes), then per frame per block float32 dx, dy, dz (i fastest).
usage:
  write_frames_bin.py out.bin --global /home/claude/frames --grid heart_grid.g --numbering node_numbering.npz \
        --csv cycle_frames.csv --period 0.857142857 --scale 1e-3
  write_frames_bin.py out.bin --blocks DIR [--period T | --no-periodic] [--scale 1]
"""
import sys, os, struct, argparse, glob, numpy as np
MAGIC = 0x4652414D

def read_g_dims(fn):
    with open(fn, 'rb') as f:
        def rec():
            n = struct.unpack('<i', f.read(4))[0]; d = f.read(n); f.read(4); return d
        nb = struct.unpack('<i', rec())[0]
        return np.frombuffer(rec(), '<i4').reshape(nb, 3).copy()

def write_header(f, dims, times, periodic, period, scale):
    f.write(struct.pack('<iiiii', MAGIC, 1, len(dims), len(times), 1 if periodic else 0))
    f.write(struct.pack('<dd', period, scale))
    f.write(np.asarray(dims, '<i4').tobytes(order='C'))      # (nblk, 3) row-major == dims(3, nblk) column-major in Fortran
    f.write(np.asarray(times, '<f8').tobytes())

def write_block(f, d):
    """d: (3, ni, nj, nk) -> dx, dy, dz each i fastest"""
    for c in range(3):
        f.write(np.ascontiguousarray(d[c].transpose(2, 1, 0), dtype='<f4').tobytes())

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('out')
    ap.add_argument('--global', dest='gdir'); ap.add_argument('--grid'); ap.add_argument('--numbering'); ap.add_argument('--csv')
    ap.add_argument('--blocks', dest='bdir')
    ap.add_argument('--period', type=float, default=0.0); ap.add_argument('--no-periodic', action='store_true')
    ap.add_argument('--scale', type=float, default=1.0)
    ap.add_argument('--zero-first', action='store_true', help='subtract frame 0 from every frame (frames relative to the first)')
    a = ap.parse_args()
    if a.gdir:
        dims = read_g_dims(a.grid)
        Z = np.load(a.numbering); inv = Z['inv']
        off = np.concatenate([[0], np.cumsum(dims.prod(1))])
        files = sorted(glob.glob(os.path.join(a.gdir, 'disp_*.npy')))
        import csv
        rows = list(csv.DictReader(open(a.csv)))
        assert len(rows) == len(files), (len(rows), len(files))
        times = [float(r['time_s']) for r in rows]
        periodic = not a.no_periodic
        with open(a.out, 'wb') as f:
            write_header(f, dims, times, periodic, a.period, a.scale)
            D0 = np.load(files[0]).astype(np.float64) if a.zero_first else None
            for n, fn in enumerate(files):
                D = np.load(fn).astype(np.float64)
                if D0 is not None: D = D - D0
                for b, (ni, nj, nk) in enumerate(dims):
                    I = inv[off[b]:off[b + 1]].reshape(ni, nj, nk)      # block node (i,j,k) -> global node
                    d = D[I]                                              # (ni, nj, nk, 3)
                    write_block(f, np.moveaxis(d, -1, 0))
                if n % 10 == 0: print(f'frame {n} t={times[n]:.5f} max|d| {np.abs(D).max():.3f}', flush=True)
        print('wrote', a.out, os.path.getsize(a.out) / 1e6, 'MB,', len(files), 'frames,', len(dims), 'blocks')
    else:
        times = np.load(os.path.join(a.bdir, 'times.npy'))
        files = sorted(glob.glob(os.path.join(a.bdir, 'frame_*.npz')))
        assert len(files) == len(times)
        z0 = np.load(files[0]); keys = sorted(k for k in z0.files if k.startswith('b'))
        dims = np.array([z0[k].shape[1:] for k in keys])
        periodic = not a.no_periodic
        with open(a.out, 'wb') as f:
            write_header(f, dims, times, periodic, a.period, a.scale)
            for fn in files:
                z = np.load(fn)
                for k in keys: write_block(f, z[k])
        print('wrote', a.out, os.path.getsize(a.out) / 1e6, 'MB,', len(files), 'frames,', len(dims), 'blocks')

if __name__ == '__main__':
    main()
