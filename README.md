# probWork3 — 심장 혈류 CFD (모핑 격자) 프로젝트

| 폴더 | 내용 |
|---|---|
| `solver/` | probCFD_부산 update 24c 솔버 소스(Fortran 90 + MPI), `Makefile`, namelist 템플릿, 문서(HANDOVER, PROJECT_RECORD), `tools/`, FSI 검증 사례 `cases/FSI_PfS2a`(샌드박스 결과는 `reference_monitor.dat`만) |
| `heart/scripts/` | 심장 격자 생성(`grid/`)·심주기(`cycle/`) 스크립트와 복원 스크립트. 절대경로 `/home/claude/hg`, `/home/claude/cyc`를 사용(`setup_resume.sh` 참고) |
| `heart/data/` | v4 최종 격자 `heart_grid.g`(143블록, 949,295셀, mm), CFDShip 형식 `.bcs`, 블록 목록, 꼭짓점 대응표, 97개 프레임 표 `cycle_frames.csv` |
| `heart/docs/` | 인계 문서(최종, v4, v1), 사용자 형식 기준 `NACA66mod_3D.bcs` |
| `tools/` | `unegg.py`(EGG 분할 압축 해제), `make_volume_frames.py`(97개 프레임 체적 모핑), `write_frames_csv.py`(프레임 STL 재생성) |

솔버 빌드: `cd solver && make clean && make` (gfortran + OpenMPI).
