# 심장 CFD 격자·심주기 작업 인계 문서 (최종)

작성일: 2026-09-28
이 문서와 `heart_cfd_resume_bundle.zip` 하나만 있으면, 새 대화방에서 지금까지의 작업 상태를 그대로 복원해 바로 이어서 작업할 수 있습니다.

---

## 0. 새 대화방에서 시작하는 방법

1. 새 대화방에 **이 문서**와 **`heart_cfd_resume_bundle.zip`**을 업로드합니다.
2. 첫 메시지 예시:

   > 첨부한 인계 문서를 먼저 읽고, 번들을 풀어서 `scripts/setup_resume.sh`로 작업 환경을 복원한 다음 `validate_resume.py`로 검증해 줘. 검증이 통과하면 ○○ 작업을 이어서 하자.

3. Claude가 할 일 (약 3분):

   ```bash
   mkdir -p /home/claude/resume && cd /home/claude/resume
   unzip -q /mnt/user-data/uploads/heart_cfd_resume_bundle.zip
   bash heart_cfd_resume_bundle/scripts/setup_resume.sh heart_cfd_resume_bundle      # 약 1분
   python3 /home/claude/hg/validate_resume.py heart_cfd_resume_bundle/data           # 약 2분
   ```

   마지막 줄에 `RESUME STATE VALID`가 나오면 복원이 끝난 것입니다.

**이 절차는 초기화된 깨끗한 환경에서 실제로 실행해 확인했습니다.** 결과는 다음과 같습니다.

| 검증 항목 | 결과 |
|---|---|
| 격자 상태 | 143블록, 949,295셀. type 92 인터페이스 798개와 물리 BC 사각형 198개가 `.bcs` 파일과 완전히 일치 |
| 유입·유출 원판 14개 | 테두리 원 오차 3×10⁻⁷ mm, 평면도 2×10⁻¹² mm |
| 경계 메시 꼭짓점 순서 | 전달된 `vertex_to_grid_nodes.txt`와 **바이트 단위로 동일** (82,012개) |
| 심실 진폭 보정 | 기록값과 차이 5.8×10⁻⁶ |
| 벽면 운동 재현 | 전달된 프레임 STL과 차이: 확장기말 0.002 mm, 수축말기 0.031 mm |
| 수축말기 morph | 뒤집힌 셀 0개, 최소 scaled Jacobian 0.0124 (원래 0.0126) |

벽면 운동의 0.03 mm 수준 차이는 **복원 오류가 아닙니다.** RBF 운동 평가가 float32로 계산되는데, float32 합산 결과는 CPU와 BLAS 커널에 따라 조금씩 달라집니다. 같은 식을 float32와 float64로 계산하면 최대 0.044 mm 차이가 나는 것을 확인했고, 관측된 차이가 그 범위 안에 있습니다. 머신과 무관하게 비트 단위로 재현해야 한다면 7.3절을 참고하세요.

---

## 1. 프로젝트 개요와 사용자 요구사항

### 1.1 목표

심장 내부 혈류(좌·우 심장과 연결된 모든 동맥·정맥) CFD 해석을 위해 다음 두 가지를 만듭니다.

1. **Multiblock PLOT3D 구조 격자**: 모든 블록이 **Block-to-Block(type 92) 경계로 연결**되고, 인터페이스의 **격자점이 정확히 일치**해야 합니다. 범위는 i, j, k 인덱스 부분범위로 지정할 수 있습니다.
2. **심주기 STL 세트**: 한 박동 주기의 벽면 형상을, CFD가 발산하지 않을 만큼 충분한 시간 단계로 나눈 것입니다. 계산 중 격자 morphing의 벽면 가이드로 씁니다.

### 1.2 사용자가 지정한 조건 (반드시 유지)

- **파일 형식**은 사용자가 준 `NACA66mod_3D.grd/.bcs`와 같아야 합니다(`docs/NACA66mod_3D.bcs` 참고).
  - `.grd`: ASCII, 첫 줄 블록 수(`%12d`), 블록별 `ni nj nk`, Fortran `0.xxxxxxxxxxxxxxxxE+yy` 형식(폭 28)으로 한 줄에 4개, i가 가장 빠르게 변함
  - `.bcs`: CRLF 줄바꿈. 블록마다 경계면 수를 적고, 각 면은 `type / 법선방향(+n=min, −n=max) / i, j, k 범위 / 0`, type 92이면 이어서 `donor 블록 / donor 법선 / donor 범위`. donor 범위의 start>end는 역방향을 뜻함
- **셀 수는 100만 개 이하**입니다(v3에서 요청). 현재 949,295개입니다.
- **격자 품질 조건** (v3에서 요청):
  1. 인접한 같은 방향 셀 길이가 비슷할 것
  2. 종횡비가 작을 것(정육면체에 가깝게)
  3. skewness가 작을 것(직교 격자에 가깝게)
- **BC 코드**는 예시를 따랐습니다: 벽 20, 유입 9, 유출 12. 사용자에게 확인 요청은 해두었지만 아직 답을 받지 못했습니다.
- **일정은 급하지 않고, 품질이 우선**이라고 사용자가 명시했습니다. 뒤집힌 셀 0개를 위해 근본 원인(토폴로지)까지 고치는 쪽을 선택했습니다.

---

## 2. 작업 이력

| 버전 | 요청 | 결과 | 핵심 |
|---|---|---|---|
| v1 (이전 방) | 해부학 기반 심장 형상과 격자 | 폴라 블록 22개, 약 152만 셀 | **블록끼리 연결되지 않음.** 혈관 블록이 방 블록 벽을 가로질러 겹쳐 있기만 함 |
| v2 | 모든 블록을 type 92 점 일치로 연결, 셀 수 유지 | 140블록, 1,451,020셀, 인터페이스 774개 | 코어 골격(박스·스윕 관)과 전역 O-grid 벽면층으로 새로 설계. 판막 사이 격벽과 IVC–삼첨판 능선을 추가 |
| v3 | 품질 조건 3가지와 100만 셀 이하 | 140블록, 944,454셀 | 벽면 방향 균일 분포, Hermite 곡선 격자선을 호길이로 분할, squircle 코어, 코어 평활화 |
| v4 | 심주기 morphing 가이드 STL | 143블록, 949,295셀, 97개 프레임 | 분지 갈림목 토폴로지 변경(벽 띠), 다중 상태 격자 최적화, 26개 위상 morph에서 뒤집힘 0개 |

v4에서 사용자가 선택한 경로:
- 격자 보강(1번안) 선택
- 격벽 여유 거리 재생성(1번안) 선택
- 격벽을 두껍고 둥글게(1번안) 선택
- "토폴로지 변경" 선택 → 이 방식으로 최종 해결

---

## 3. 현재 결과물 (v4)

### 3.1 격자 (`data/heart_grid.g`, 전체 결과물에는 `.grd`도 있음)

| 항목 | 값 |
|---|---|
| 블록 / 셀 | 143 (좌심 93, 우심 50) / 949,295 |
| type 92 인터페이스 | 798개, 격자점 불일치 0.000 mm, 모두 same-order, transpose 없음 |
| scaled Jacobian | 최소 0.046, p1 0.456, 중앙값 0.989 |
| 면 비직교각 | 중앙값 5.9°, p95 33.4° |
| 종횡비 | 중앙값 1.79, p99 6.5 (최댓값 422: 분지 입구·격벽 부근의 셀 35개) |
| 인접 셀 크기비 | p95 1.10, p99 1.31 |
| 유체 체적 (휴지 형상) | 좌 428 mL, 우 400 mL |
| 벽면 첫 셀 | 균일 분포. 대동맥 기준 y+ 약 8–9 (표준 벽함수에는 부적합) |

블록 번호와 이름은 `data/block_list.txt`에 있습니다. 이름 형식은 코어 블록이 `LV`, `AO`, `pv_left_inferior` 등이고, 벽면층 블록이 `<코어>.<면>.bl<번호>`입니다.

### 3.2 심주기

- **생리 설정**: 70 bpm(주기 857 ms). LV 133→50 mL, RV 140→62 mL로 박출량 −83 / −78 mL이며 정확히 보정. 방실판 하강 약 12 mm, 비틀림 ±6°, 심방 등방 변화, 대혈관 맥압 팽창.
- **프레임**: 97개, 간격 1.2–9.8 ms.
  - 선정 기준: 인접 프레임을 선형 보간했을 때 오차 ≤ 벽면 첫 셀 높이의 20%(하한 0.15 mm), 간격 ≤ 10 ms
  - 더 촘촘하게 하려면 기준을 10%로 낮추면 됩니다(271개 프레임).
- **벽면 검증 (전 프레임)**: 삼각형 뒤집힘 0개, 최소 면적비 0.24, 좌·우 최소 간격 0.65 mm, 최대 변위 16.0 mm.
- **격자 morph 검증**: 확산 morphing(벽 거리^−1 강성) 기준, 26개 위상에서 뒤집힘 0개. 최솟값은 0.011(s=0.36)입니다.
- **유입·유출 면**: 모두 평면 유지. 관상동맥 출구 2개만 평행이동하고 나머지는 고정입니다.
- **꼭짓점 대응**: 모든 프레임의 꼭짓점·삼각형 순서가 같습니다. `vertex_to_grid_nodes.txt`가 꼭짓점(0부터) → `블록:i,j,k`(1부터)를 알려줍니다.

### 3.3 이전 방의 전체 결과물 (`/mnt/user-data/outputs`, 이 번들에는 일부만 포함)

| 파일 | 설명 | 번들 포함 |
|---|---|---|
| `heart_cycle_stl_97frames.zip` (346 MB) | 97개 프레임 전체 | 없음. 검증용으로 프레임 000·031의 꼭짓점 좌표만 포함. 전체는 약 5분에 재생성 가능(6.3절) |
| `heart_grid.grd` (100 MB) | ASCII 격자 | 없음. `export.py`로 재생성(6.4절) |
| `heart_grid.g`, `.bcs`, `block_list.txt`, `vertex_to_grid_nodes.txt`, `cycle_frames.csv` | — | 포함 |
| `heart_reference_stl_v4.zip` | 휴지 형상 STL | 없음. `export.py`로 재생성 |
| `HANDOFF_심장CFD_v4_심주기.md` | v4 기술 문서 | `docs/`에 포함 |

---

## 4. 번들 구성 (`heart_cfd_resume_bundle.zip`)

```
heart_cfd_resume_bundle/
├─ HANDOFF_심장CFD_작업인계_최종.md     ← 이 문서
├─ scripts/
│  ├─ setup_resume.sh      환경 복원 (패키지 설치 → 스크립트 배치 → 상태 복원 → 경계 메시 → 보정값)
│  ├─ restore_state.py     heart_grid.g + .bcs + block_list → assembled_s5.pkl (자체 검증 포함)
│  ├─ validate_resume.py   원래 세션과 같은지 검증 (보정값, 프레임 재현, 수축말기 morph)
│  ├─ verify_bcs.py        .grd/.bcs 파일만 읽어 점 일치와 면 지정을 독립 검증
│  ├─ grid/   → /home/claude/hg 로 복사됨
│  └─ cycle/  → /home/claude/cyc 로 복사됨
├─ data/
│  ├─ heart_grid.g  heart_grid.bcs  block_list.txt  vertex_to_grid_nodes.txt
│  ├─ cycle_frames.csv  ES_FRAME.txt  reference_frames.npz (프레임 000과 031의 꼭짓점 좌표, 검증용)
│  └─ morph_validation_phases.json  quality_reference_grid.txt
└─ docs/
   ├─ HANDOFF_심장CFD_v4_심주기.md       v4 기술 문서
   ├─ HANDOFF_심장CFD모델_v1.md          v1 원본 인계 문서 (해부학 파라미터 설명)
   └─ NACA66mod_3D.bcs                   사용자가 준 형식 기준 파일 (.grd는 266 MB라 제외)
```

**경로 주의**: 스크립트는 `/home/claude/hg`(격자)와 `/home/claude/cyc`(심주기)를 절대경로로 사용합니다. `setup_resume.sh`가 이 구조를 만들어 줍니다.

---

## 5. 코드 구조

### 5.1 격자 (`/home/claude/hg`)

| 파일 | 역할 |
|---|---|
| `params.py` | v1 해부학 파라미터(좌표 mm, 방 중심과 반경, 혈관 중심선, 판막 위치) |
| `paths.py` | 중심선(Hermite), implicit 벽(smooth union, 격벽 `cuts`, 좌·우 상호 배제) |
| `core.py` | 코어 블록 생성: 스윕 관 `tube`(squircle 단면 `SQ=0.30`, 단면 Laplacian 평활화, 끝단 고정), `sbox`(구형화 박스), 면 법선, 부분범위 선택 |
| `heart.py` | 심장 코어 골격 설계. `RES` 파라미터, 좌/우 조립, 분지 부착, 대동맥 비틀림, 격벽 `curtain`, 벽 띠(`Ngap`) |
| `oblayer.py` | 전역 O-grid 벽면층: 방향장 평활화, 벽까지 광선, 두께 보정, 호길이 Hermite 기둥(`kk0=0.35`, 뒤집히는 절점만 곡률 완화) |
| `mb.py` | 점 일치 인터페이스 자동 탐지(`connectivity`), BC 분해(`bc_surfaces`), `.bcs`/`.grd`/`.g` 출력, Jacobian, 비직교각 |
| `final_dom.py` | 도메인별 생성. 인자: `L|R 층수 축간격 좌우배제간격 해상도배율` → `final_L.pkl`/`final_R.pkl` |
| `assemble.py` | 조립, 유입·유출 자동 태그, 전 항목 검증 → `assembled.pkl` |
| `optimize_grid.py` | 정적 품질 최적화(내부 격자점, 경계 고정) |
| `optimize_motion.py` | **다중 상태 최적화**: 기준 격자와 여러 morph 상태의 scaled Jacobian을 동시에 개선. 벽 격자점은 원래 벽면 삼각형 위에서만 미끄러짐. 유입·유출은 고정. 기준 격자 품질 보호(`REFMIN`) |
| `numbering.py` | 전역 격자점 번호 고정. 기준 파일 `assembled_base.pkl`, 캐시 `node_numbering.npz` |
| `export.py` | `GRID=... python3 export.py OUTDIR` → `.grd`/`.g`/`.bcs`/`block_list`/STL |
| `quality.py`, `qstat.py` | 좌·우 간섭과 품질 지표 |

### 5.2 심주기 (`/home/claude/cyc`)

| 파일 | 역할 |
|---|---|
| `deform_v1.py`, `geom.py` | v1의 심주기 곡선(LV/RV/LA/RA 용적, 대혈관 팽창, 위상 이름)과 가중 함수 |
| `cycle.py` | 해석적 방 운동 `Motion`(도메인별 분리, 관상동맥 출구 평행이동, 비틀림 `twist_deg`, 장축 계수 `k_lv`) |
| `bmotion2.py` | 운동 지정 영역 선택(심실 몸체, 심방 위쪽, 60 mm 밖), 갈림목 주변 해제(`crotch_r=14`), 좌·우 접촉 결합, 맥압 팽창 평활화 |
| `bmotion3.py` | **최종 벽면 운동**: 자유 벽을 3D RBF(r³+선형, 제어점 좌 2600 / 우 1600)로 보간, 잔차는 벽면 조화 확장으로 보정, 원거리 감쇠 25–60 mm |
| `crotch.py` | 갈림목 좌표: LV·RV 기저부 벽 띠 중앙선, PA 끝면, 삼첨판 하부(IVC 능선) |
| `bmesh.py` | 격자 경계 → 삼각형 메시(`bmesh.npz`)와 꼭짓점 대응표 |
| `diffmorph.py` | 확산 morphing 검증기(대수 다중격자 `pyamg`, 프레임당 약 50초) |
| `outer1.sh <grid이름> <태그>` | 격자가 바뀐 뒤 경계 메시 재생성과 보정 |
| `save_morph.py <grid.pkl> <out.npz> <위상들>` | 위상별 morph 상태를 `out_sXXXX.npy`로 **즉시** 저장 (optimize_motion 입력) |
| `fullcycle.py`, `check_morph.py` | 위상별 morph 검증 (`check_morph`는 뒤집힌 셀이 있는 블록 이름까지 표시) |
| `sample_cycle.py`, `pick_frames.py`, `write_frames.py` | 700개 표본 → 적응형 프레임 선정 → STL과 CSV 출력 |

### 5.3 주요 파라미터 (최종값)

- **`heart.RES`**: `lean=0.7`(LV 기저부 승모판·LVOT 약 35° 벌림), `lean_rv=0.4`, `Ngap=3`, `Ngap_rv=3`, `Ngap_pa=2`(벽 띠), `cclear=3.0`, `ctau=1.25`, `ck=1.5`(격벽: 코어 여유 3 mm, 두께 약 2.5 mm), `isthmus=11`, `vein_tilt=20`, `kblend=6`, `beta_a=0.62`
- **생성 명령**: `final_dom.py L 9 1.0 0.45 0.85` 와 `R 9 1.0 0.45 0.85` (벽면층 9층, 축 간격 1.0 mm, 좌·우 배제 0.45 mm, 해상도 배율 0.85)
- **벽면 운동**: `far_mm=60`, `crotch_r=14`, `fade_far=(25,60)`, `nctrl=(2600,1600)`, `wexp=2`
- **보정값**: `a_lv=1.0210203528`, `a_rv=0.9938990474`
- **morph 검증**: `DiffMorph(alpha=1.0)`
- **프레임 선정**: `beta=0.2`, `h_floor=0.15 mm`, `dt_max=10 ms`

---

## 6. 자주 하는 작업의 절차

### 6.1 격자 상태 검사

```bash
cd /home/claude/cyc && GRID=/home/claude/hg/assembled_s5.pkl python3 check_morph.py 0.2 0.33 0.36 0.45
```

### 6.2 격자를 수정했을 때 (다중 상태 최적화 외부 반복 1회)

```bash
cd /home/claude/cyc
bash outer1.sh <새grid> <태그>                                            # 경계 메시와 보정 갱신
python3 save_morph.py /home/claude/hg/<새grid>.pkl st_<태그>.npz 0.2 0.33 0.45   # 위상당 약 50초
cd /home/claude/hg
python3 optimize_motion.py <새grid>.pkl <다음grid>.pkl 0.35 100 230 /home/claude/cyc/st_<태그>_s0200.npy ...
```

- 수렴하면 수축기 0.30–0.42 구간을 0.01–0.02 간격으로, 전체 주기를 0.05 간격으로 검증합니다.
- **블록 구조(블록 수·차원)가 바뀌면** 먼저 `assembled_base.pkl`을 새 격자로 바꾸고 `node_numbering.npz`를 지워야 합니다.

### 6.3 심주기 프레임 재생성

```bash
cd /home/claude/cyc
python3 sample_cycle.py 700 0 700      # 약 160초 → Dsamples.npy
python3 pick_frames.py                 # 프레임 후보 → frame_sets.json
python3 write_frames.py                # 출력 경로(write_frames.py 안) 확인 후 실행
```

### 6.4 격자 파일 출력과 독립 검증

```bash
cd /home/claude/hg && GRID=assembled_s5.pkl python3 export.py /home/claude/out
python3 verify_bcs.py /home/claude/out/heart_grid.grd /home/claude/out/heart_grid.bcs
```

---

## 7. 교훈과 주의 사항 (반드시 읽을 것)

### 7.1 실패했던 접근 (반복하지 말 것)

- **분지부의 칼날형 분할이 수축기 격자 뒤집힘의 근본 원인이었습니다.** 두 관이 부모 면의 한 격자선을 공유하는 구조를 말합니다. 다음 조치는 모두 **효과가 없었습니다**: 격벽 여유 거리·두께·둥글림 조정, 코어 각도 변경, 갈림목 운동을 보간이나 강체로 바꾸기, 비틀림 제거, 격벽 양면 결합, 방실판 하강 완화. **부분범위 사이에 2–3셀 벽 띠를 두니 해결됐습니다.**
- **벽면 위에서만 매끄러운 운동**(벽면 조화 확장)은 벽면 품질은 좋지만, 공간적으로 가까운 두 벽이 서로 다르게 움직여 내부 격자를 망가뜨립니다. **반드시 3D 보간(RBF)**을 써야 합니다.
- **갈림목 주변만 강체 운동**시키면 수축하는 기저부와 어긋나 322셀이 뒤집혔습니다.
- **눌린 벽면층 기둥을 벽 쪽으로 두껍게** 하면 기준 격자에 새 뒤집힘이 생깁니다.
- **RBF morph 검사는 신뢰할 수 없습니다.** 실제 CFD morpher와 같은 **확산 morphing으로 검증**해야 합니다.
- 코어 평활화 때 **다른 블록과 맞닿는 첫·끝 단면까지 움직이면** 점 일치가 깨집니다. 끝단은 고정하고 부근은 점진적으로 평활화해야 합니다.

### 7.2 도구 환경 제약

- 명령 한 번의 실행 한도는 **300초**입니다. 긴 작업은 나누고 결과를 즉시 파일로 저장해야 합니다(`save_morph.py`가 위상별 `.npy`를 바로 쓰는 이유).
- **백그라운드 작업은 응답이 끝나면 종료**됩니다. 한 호출 안에서 끝나도록 해야 합니다.
- 메모리는 약 **3 GB**입니다. 좌·우 도메인을 동시에 생성하면 메모리 부족으로 강제 종료됩니다.
- 격자가 바뀌면 벽면 메시(`bmesh.npz`)도 **반드시** 다시 만들어야 합니다(`outer1.sh`). 그렇지 않으면 `DiffMorph`에서 `AssertionError`가 납니다.

### 7.3 재현성

- 격자, 경계 메시, 꼭짓점 순서는 파일에서 **완전히** 복원됩니다.
- 벽면 운동은 float32 RBF 평가 때문에 **머신에 따라 최대 약 0.05 mm** 다를 수 있습니다. 격자 유효성에는 영향이 없음을 확인했습니다.
- 머신과 무관하게 재현해야 하면 `bmotion3.py`에서 `E`를 `float64`로 저장하도록 바꾸면 됩니다. 메모리가 약 1 GB 늘고, 그 뒤로 만든 프레임은 전달된 프레임과 최대 0.044 mm 다릅니다.
- `final_dom.py`부터 다시 생성하면 v4의 **최적화 전** 격자가 나옵니다. v4 최종 격자(`assembled_s5`)는 외부 반복 최적화를 여러 번 거친 결과이므로, **처음부터 재생성하지 말고 파일에서 복원**하세요.

---

## 8. 미해결 사항과 다음 작업 후보

1. **BC 코드 확인**: 벽 20, 유입 9, 유출 12가 사용자 솔버에 맞는지 아직 답을 받지 못했습니다. 바꾸려면 `mb.py`의 `BCCODE`만 수정하고 `export.py`로 다시 출력합니다.
2. **y+**: 현재 약 8–9입니다. 벽함수를 쓰면 `final_dom.py`의 6번째 인자(벽면 집중계수)나 층수를 조정해야 하는데, **격자가 바뀌므로 심주기 최적화부터 다시** 해야 합니다.
3. **관상동맥 출구가 움직임**: 솔버가 출구 이동을 허용하지 않으면 짧은 연장관을 붙이는 방안이 있습니다.
4. **솔버 연동 시험**: 사용자 솔버의 morpher가 확산 방식이 아니면 결과가 다를 수 있습니다. 문제가 생기면 해당 위상과 위치를 받아 `optimize_motion.py`에 그 위상을 추가해 최적화합니다.
5. **판막 소엽, 변형 시계열의 다른 심박수, 병변 모델** 등 v1 인계 문서의 미해결 항목은 아직 다루지 않았습니다.
6. **더 촘촘한 프레임 세트**가 필요하면 `pick_frames.py`에서 `beta=0.1`로 설정합니다(271개 프레임).
