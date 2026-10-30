#!/bin/bash
# 새 대화방에서 심장 CFD 작업 환경을 복원한다.
# 사용법:  bash setup_resume.sh <resume 번들을 푼 폴더>
#   - /home/claude/hg  (격자 스크립트 + 복원된 격자 상태)
#   - /home/claude/cyc (심주기 스크립트 + 경계 메시 + 보정값)
set -e
SRC=$(cd "${1:-.}" && pwd)
echo "== 1/5 패키지 설치"
pip install -q --break-system-packages numpy scipy trimesh pyamg matplotlib manifold3d 2>&1 | tail -1 || true
python3 -c "import numpy,scipy,trimesh,pyamg; print('   numpy',numpy.__version__,'scipy',scipy.__version__,'trimesh',trimesh.__version__,'pyamg',pyamg.__version__)"
echo "== 2/5 스크립트 배치"
mkdir -p /home/claude/hg /home/claude/cyc
cp "$SRC"/scripts/grid/*.py /home/claude/hg/
cp "$SRC"/scripts/cycle/*.py "$SRC"/scripts/cycle/*.sh /home/claude/cyc/
cp "$SRC"/scripts/restore_state.py "$SRC"/scripts/validate_resume.py "$SRC"/scripts/verify_bcs.py /home/claude/hg/
echo "== 3/5 격자 상태 복원 (heart_grid.g + .bcs + block_list.txt -> assembled_s5.pkl)"
python3 /home/claude/hg/restore_state.py "$SRC/data"
echo "== 4/5 경계 메시 + 꼭짓점 대응표 재생성"
cd /home/claude/cyc
GRID=/home/claude/hg/assembled_s5.pkl OUT=/home/claude/cyc/bmesh.npz MAP=/home/claude/cyc/vertex_to_grid_nodes.txt python3 bmesh.py
if cmp -s /home/claude/cyc/vertex_to_grid_nodes.txt "$SRC/data/vertex_to_grid_nodes.txt"; then
  echo "   꼭짓점 순서: 전달본과 동일 (82,012개)"
else
  echo "   !! 꼭짓점 순서가 전달본과 다름 - 확인 필요"; exit 1
fi
echo "== 5/5 심실 진폭 보정값 기록"
python3 - <<PY
import pickle
cal={'a_lv': 1.0210203528404236, 'a_rv': 0.9938990473747253}   # 원래 세션의 최종 보정값 (격자 s5)
for f in ('calib_final.pkl','calib_r.pkl'): pickle.dump(cal,open('/home/claude/cyc/'+f,'wb'))
print('   calib_final.pkl / calib_r.pkl :',cal)
PY
cp "$SRC/data/cycle_frames.csv" /home/claude/cyc/ 2>/dev/null || true
echo "복원 완료. 검증: python3 /home/claude/hg/validate_resume.py $SRC/data"
