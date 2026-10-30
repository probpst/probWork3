# 심장 혈류 CFD 모델링 — 작업 인계 문서 (HANDOFF)

> 이 문서 하나만 있으면 다른 세션/다른 사람이 바로 이어서 작업할 수 있도록 정리한 것입니다.
> 새 대화방에서는 **이 파일과 `scripts/` 폴더(9개 .py)만 업로드**하면 전체 재현이 가능합니다.

---

## 0. 한 줄 요약

175 cm 표준 성인의 **심장 + 대혈관 혈액 유체영역(blood pool)** 을 파라메트릭하게 생성하고,
여기에 **판막 지오메트리 / PLOT3D 정렬격자 / 한 심주기 변형 STL 시계열** 까지 만든 작업.
모두 Python(numpy·scipy·trimesh·manifold3d)으로 절차적 생성했으며,
환자 영상 세그멘테이션이 아닌 **이상화(idealized) 모델** 입니다.

---

## 1. 최종 산출물 목록

작업 디렉터리 루트 = `/mnt/user-data/outputs/`

| 경로 | 크기 | 내용 |
|---|---|---|
| `heart_cfd_package.zip` | 148 MB | 아래 `heart_cfd/` 전체 |
| `heart_grid_plot3d.zip` | 51 MB | PLOT3D 격자만 따로 |
| `heart_cycle_stl.zip` | 125 MB | 심주기 변형 STL 40프레임 |
| `heart_3d_render.png` | — | 전체 형상 3면 렌더 |

### 1.1 `heart_cfd/` — 정지 형상 (t/T = 0, 이완말기)

```
heart_cfd/
├── README.md                         ← 상세 문서 (8개 절)
├── model_metadata.json               ← 모든 패치 면적/중심/법선/용적
│
├── left_heart_fluid.stl              346,172 tri  424.9 mL  watertight
├── right_heart_fluid.stl             300,052 tri  391.7 mL  watertight
├── heart_full_merged.stl             646,224 tri  816.6 mL  ★통합 STL
├── heart_full_merged_with_valves.stl 859,140 tri   (+개방 첨판 포함)
│
├── left_heart/    inlet_pv_*.stl ×4, outlet_*.stl ×6, wall.stl
├── right_heart/   inlet_svc/ivc.stl, outlet_rpa/lpa.stl, wall.stl
├── valve_planes/  aortic / mitral / pulmonary / tricuspid _annulus_plane.stl
├── valves/        {aortic,pulmonary,mitral,tricuspid}_{open,closed}/ 첨판 solid
│
├── left_heart_systole.stl / left_heart_diastole.stl
├── right_heart_systole.stl / right_heart_diastole.stl   ← 개방 첨판 boolean 제거본
│
├── grid/          heart_grid.g (PLOT3D binary), heart_grid.fmt (ASCII),
│                  boundary_conditions.txt, grid_report.json, grid_preview.png
│
├── preview_overview.png / preview_valves_axial.png / preview_valves_side.png
└── scripts/       geom.py params.py build.py assemble.py valves.py
                   make_valves.py export_all.py grid.py make_grid.py
```

### 1.2 `heart_cycle/` — 심주기 변형 시계열

```
heart_cycle/
├── README_CYCLE.md
├── heart_cycle_000.stl … heart_cycle_039.stl   (40 프레임, 각 147,596 tri)
├── base_reference_mesh.stl                     (변형 기준 메시)
├── cycle_metadata.json                         (프레임별 시간/위상/체적/최대변위)
├── cycle_curves.png / cycle_phases.png
└── scripts/  deform.py params.py geom.py
```

---

## 2. 좌표계 · 단위 (전 산출물 공통)

| 축 | 방향 |
|---|---|
| **+x** | 환자 **좌측** |
| **+y** | **전방** (흉골 쪽) |
| **+z** | **상방 = 바닥 기준 높이** |

- 단위 **mm**, 원점 = 정중시상면 × 중간관상면 근사 × 바닥
- 주요 높이: 심첨 z≈1248 / 승모판륜 1322 / 대동맥판륜 1332 / 폐동맥판륜 1356 /
  대동맥궁 정점 1421 / 하행대동맥 출구 1180
- 전체 bounding box: x ∈ [−62.6, 80.4], y ∈ [−69.2, 70.8], z ∈ [1180, 1480]

---

## 3. 코드 구조 (scripts/)

의존성: `pip install trimesh manifold3d shapely networkx numpy scipy rtree triangle`

| 파일 | 역할 |
|---|---|
| `params.py` | ★ **모든 해부학적 치수가 여기 한 곳에.** 중심선 제어점·반경·타원체 파라미터 |
| `geom.py` | 기하 커널: Catmull-Rom 리샘플, parallel-transport 프레임, 스윕/로프트, Delaunay 평면 캡, manifold3d 불리언 |
| `build.py` | 구성요소별 메시 생성 + **간섭 검사(`FORBIDDEN` 리스트)** |
| `assemble.py` | 불리언 합집합 + 경계 패치 분리 |
| `valves.py` | 판막 첨판 파라메트릭 생성기 (반월판 / 방실판) |
| `make_valves.py` | 4개 판막 × 개폐 2상태 구동 |
| `export_all.py` | 전체 파이프라인 (유체영역·패치·판륜면·첨판·위상별 도메인) |
| `grid.py` | 정렬격자 커널 (폴라 블록, 회전타원체 블록, RV 초승달 H블록, y+ 기반 반경분포) |
| `make_grid.py` | 22블록 조립 + 품질검사 + PLOT3D 쓰기 |
| `deform.py` | 심주기 변형장 + 시계열 STL 생성 |

### 재실행 명령

```bash
# 1) 정지 형상 전체 (고해상도, 약 3~4분)
python3 export_all.py highres /출력/heart_cfd 121

# 2) 저해상도 + 간섭 검사만 (치수 수정 후 반드시 먼저 실행)
python3 build.py lowres          # -> "clearance violations: none" 확인

# 3) 정렬격자
python3 make_grid.py /출력/heart_cfd/grid

# 4) 심주기 변형 시계열
python3 deform.py <기준.stl> <출력폴더> <프레임수> <bpm>
```

---

## 4. 형상 모델 — 설계와 검증

### 4.1 구성

- **좌심계**: 폐정맥 4개 → 좌심방 → 승모판륜 → 좌심실 → LVOT → 대동맥근(Valsalva 동 3엽)
  → 상행 → 궁 + 3분지(완두동맥간·좌총경·좌쇄골하) → 하행흉부대동맥, + 관상동맥 개구부 2개
- **우심계**: SVC/IVC → 우심방 → 삼첨판륜 → 우심실(중격을 감싼 **초승달 단면**)
  → RVOT → 주폐동맥 → LPA/RPA
- 좌/우심계는 **분리된 2개 solid**. 교차 체적 0.0000 mL 검증 완료.

### 4.2 생성 방식

- 혈관·판막트랙 = 중심선 스플라인 + 반경 프로파일 **스윕**
- 좌심실 = 절단된 **회전타원체 로프트** (자오선 호길이 균등 샘플링)
- 우심실 = (RV 외곽 회전타원체) − (LV 심외막 회전타원체) **불리언 차집합** → 초승달 단면
- 심방 = 아이코스피어 기반 타원체
- 대동맥근 = 3엽(clover) 단면 프로파일 `r(θ) = R(1 + k cos3(θ−0.35))`, k: 판륜 0 → 동 0.20 → STJ 0

### 4.3 검증값 (전부 정상 성인 범위 안)

| 항목 | 모델 | 정상 |
|---|---|---|
| LV EDV | 133 mL | 106–160 |
| RV EDV | 140 mL | 120–180 |
| LA / RA | 46 / 52 mL | 40–60 |
| 대동맥판륜 | Ø24.9 mm (4.88 cm²) | 23–25 |
| 승모판륜 | 8.04 cm² (Ø32.0) | 8–10 |
| 폐동맥판륜 | Ø23.0 mm | 22–24 |
| 삼첨판륜 | Ø35.0 mm (9.61 cm²) | 33–40 |
| 상행대동맥 | Ø30 mm | 28–34 |

메시: 평균 edge 0.68–0.77 mm, 최대 1.2 mm, watertight, 단일 body.

### 4.4 ⚠ 간섭 회피 — 가장 손이 많이 간 부분

심장은 실제로 구조물끼리 맞닿아 있어서, 치수를 조금만 바꿔도 겹칩니다.
`build.py`의 `FORBIDDEN` 리스트가 **서로 닿으면 안 되는 쌍**을 정의하고
불리언 교차 체적으로 검사합니다. 실제로 해결에 반복이 필요했던 쌍:

- 상행대동맥 ↔ 주폐동맥 (서로 붙어 지나감) → 대동맥을 우측으로, MPA를 좌전방으로
- 상행대동맥 ↔ SVC → SVC를 x = −44로 이동
- 우폐동맥(RPA) ↔ 상행대동맥 (RPA가 대동맥 뒤를 지남) → RPA를 y = −30까지 후방으로
- 좌심방 ↔ 하행대동맥 (실제로 3~5 mm 간격) → 미세 조정
- 대동맥근 ↔ 우심방 / 삼첨판 트랙 (막성중격 부위)
- 좌관상동맥 개구부 ↔ RVOT
- 우폐정맥 ↔ 하행대동맥 → 폐정맥을 전방 경로로 재라우팅

**치수 수정 시 반드시 `python3 build.py lowres` 로 `violations: none` 확인할 것.**

### 4.5 판막

- 반월판(대동맥·폐동맥): 첨판 3개. 부착선은 scalloped(`z = h_c·s²`), 교련부 높이 15/13 mm,
  폐쇄 시 중심 Y형 접합, 개방 시 동벽에 밀착. 두께 0.5–0.6 mm
- 방실판(승모 2엽 / 삼첨 3엽): saddle 판륜, 폐쇄 시 접합선 tenting, 개방 시 깔때기형
- **부착선을 루멘 벽보다 1.2 mm 바깥으로** 빼서 불리언 시 누설 틈이 없게 함
- 첨판 부착 반경은 레이캐스팅이 아니라 **해석적 대동맥근 프로파일** 사용
  (레이캐스팅은 LVOT/관상동맥 개구부로 빠져나가 첨판이 튀는 버그가 있었음)

---

## 5. PLOT3D 정렬격자

**22 블록 / 1,518,528 셀 / 음의 체적 셀 0개 / 단위 mm**

- `heart_grid.g` : Fortran unformatted, record marker 있음, little-endian, double, 3D whole, IBLANK 없음
- `heart_grid.fmt` : 동일 데이터 ASCII (검증용)
- ParaView: PLOT3D Reader → Binary On, Has Byte Count On, Multi Grid On, IBlanking Off, Little Endian

### 5.1 토폴로지

- **폴라 블록 19개** (혈관·판막트랙·대혈관): i=중심선, j=반경(j=1이 축 → **POLE BC**),
  k=원주(k=1과 k=NK 일치 → 주기 1-to-1), j=NJ가 벽
- **좌/우심방·좌심실**: 해석적 내벽 폴라 블록, 극점/심첨은 `SINGULAR_POINT`
- **우심실(블록 17)**: 초승달 단면이라 폴라 불가 → **H형 블록**
  (i=기저→심첨, j=중격벽→자유벽 **양쪽 클러스터링**, k=초승달 둘레)

### 5.2 블록 연결

| 종류 | 개수 |
|---|---|
| `MATCHED 1-to-1` | **1** (블록 04↔05 대동맥 협부) |
| `PERIODIC` (자기 연결) | 19 |
| `PATCHED` (기하 정합, 인덱스 불일치 → 보간 필요) | **15** |

> ⚠ 심장은 Y자 분지 구조라 폴라 블록만으로는 접합부를 점 일치로 만들 수 없습니다.
> 완전 1-to-1 토폴로지가 필요하면 각 단면을 butterfly(O-H) 5블록으로 쪼개고
> 분지부에 전용 접합 블록을 넣어야 하며, 이는 ICEM CFD/Pointwise 대화형 블로킹 영역입니다.

### 5.3 ★ y+ 20–30 과 150만 셀의 충돌 — 해결 방법 (핵심 설계 결정)

혈액 ν = 3.30e-6 m²/s (물의 3.3배), Re 5,000–11,000 으로 **낮음**.
`u_τ = U√(C_f/2)`, `C_f = 0.079 Re^-0.25` 로 계산하면 y+=25 가 요구하는
벽면 첫 셀 높이는 **2~9 mm** — 혈관 반경과 맞먹습니다.

통상적인 "벽에서 멀어질수록 커지는" 분포를 쓰면 반경 방향에 5~6셀밖에 안 들어가
150만 셀을 채울 수 없습니다.

**해결**: 반경 분포를 **벽면 셀이 가장 크고 축 쪽으로 수축하는 형태(ratio ≈ 0.85)** 로 구성.
→ 벽면 첫 셀은 y+ 조건을 정확히 만족하면서 반경 방향 7~14 격자점 유지.
인접 셀 크기비 1.16–1.29 로 완만. (`grid.py`의 `_geom_ratio`가 g<1 을 허용)

달성 y+: 대동맥 25.3 / 하행대동맥 24.9 / 주폐동맥 25.5 / 승모·삼첨 트랙 24.7–25.0 /
좌심실 22.2 / 심방 18 / 폐정맥·대정맥·관상동맥 5–14 (Re 340–2,100 층류라 물리적으로 불가능)

→ **k-ω SST 또는 k-ε + wall function 전제.** low-Re 모델이면 `grid.py`의 `YPLUS=1` 로 재생성.

### 5.4 기타 격자 이슈

- 대동맥궁처럼 곡률 큰 구간은 `r ≤ 0.82/κ` 로 반경 제한 → 오목면 셀 반전 방지
  (궁 안쪽에서 수 mm 의 STL 대비 형상 오차 발생)
- RV 초승달 끝(horn)은 두께 0 이 되므로 θ 범위를 3.5% 줄이고 축방향 끝을 3.5% 잘라냄
- 22 블록 체적 합 762 mL = STL 816.6 mL 의 93% (차이 = 접합부 중첩 제거분, RV horn/심첨, 심방이)

---

## 6. 심주기 변형 시계열

**70 bpm, T = 0.857143 s, 40 프레임, Δt = 21.43 ms, 각 147,596 tri**

### 6.1 핵심 설계

- **하나의 기준 메시 절점을 해석적 변위장으로 움직임** (프레임마다 불리언을 다시 돌리지 않음)
  → 모든 프레임의 삼각형 순서·꼭짓점 순서가 동일 → 두 프레임 빼면 그대로 변위장
  → trimesh는 `process=False` 로 읽을 것 (기본값은 절점 병합 순서가 달라질 수 있음)
- 가중장이 전부 공간상 매끄러운 함수(smoothstep) → 구조 경계에서 메시가 찢어지지 않음

### 6.2 운동 모델

| 구조 | 운동 |
|---|---|
| 좌/우심실 | 심첨 고정 **장축 단축 15%** + **반경 수축 33%** + **비틀림 12°** |
| 심방 | 중심 기준 등방 스케일링 |
| 대동맥·폐동맥 | 중심선 기준 반경 맥압 팽창 (×1.062 / ×1.050) |
| 폐정맥·대정맥·궁분지 | 거의 강체 |

체적: LV 133→50 (EF 62%), RV 140→62 (EF 56%), LA 26↔46, RA 29↔52. CO ≈ 5.8 L/min

위상: 0–0.05 등용성수축 / 0.05–0.36 박출 / 0.36–0.43 등용성이완 /
0.43–0.58 급속충만 / 0.58–0.85 diastasis / 0.85–1.00 심방수축

### 6.3 체적 보정

가중장 기반 변형은 목표 체적을 그대로 주지 않으므로,
**심실 변형 진폭을 실제 메시 체적 응답으로 이분법 역보정**함
(보정계수 LV 1.019 / RV 0.993 → 목표 −83.0 / −78.0 mL 정확히 달성).
`Deformer.calibrate()` 참조. 좌/우 solid 식별은 `contains(O_LV)` 로 판정
(중심점 거리 비교는 오판했음 — 주의).

### 6.4 프레임 수 근거

프레임당 최대 절점 이동량 **2.94 mm** (급속충만기), 대부분 0.1–2.2 mm.
메시 평균 변 1.45 mm, 격자 벽면 셀 2.2–4.4 mm 이므로 프레임당 격자 1~2칸.
`cycle_metadata.json`의 `max_vertex_step_mm` 로 구간별 확인 가능.

---

## 7. CFD 적용 시 주의 (사용자에게 반복 강조한 것)

1. **입출구에 유량을 규정하지 말 것.** 벽이 움직여 체적이 변하므로 유량이 자동 생성됨.
   유량 BC를 같이 걸면 과구속 → 확실히 발산. 입구 압력 + 출구 Windkessel 권장.
2. **프레임 사이는 3차 스플라인(PCHIP/Catmull-Rom) 보간.** 선형보간은 벽면 속도가
   구간별 상수 → 프레임 경계에서 속도 불연속 → 압력 스파이크 → 발산 원인.
   주기 경계(039 다음 000) 적용 필수.
3. **판막 닫기.** 시계열에 판막이 없으므로 수축기엔 승모·삼첨판 평면을,
   이완기엔 대동맥·폐동맥판 평면을 배플로 닫지 않으면 역류 발생.
   `heart_cfd/valve_planes/*.stl` 사용.
4. 혈액 물성: ρ=1060 kg/m³, μ=3.5e-3 Pa·s. 비뉴턴이면 Carreau
   (μ∞=0.0035, μ0=0.056, λ=3.313 s, n=0.3568)
5. 유량 분배: 폐정맥 RSPV 27/RIPV 25/LSPV 25/LIPV 23 %,
   대동맥 하행 60/완두 18/좌총경 6/좌쇄골하 11/관상 4–5 %,
   SVC 35 : IVC 65, RPA 55 : LPA 45
6. 상행대동맥 peak Re 4,000–5,000, Womersley α 15–18 → **난류 모델 필요**

---

## 8. 모델의 한계 (매번 명시한 내용)

**적합**: CFD 방법론 개발·검증, 격자 수렴성/난류모델 비교, 판막 설계 초기 스크리닝, 교육·벤치마크

**부적합 / 누락**
- ❌ **임상 의사결정, 환자 특이적 치료 계획 — 절대 사용 금지**
- 유두근·건삭(chordae)·심실 육주(trabeculae) 없음 → LV 내 와류 구조와 WSS 분포가 실제와 다름
- **좌심방이(LAA) 없음** → 심방세동 혈전·색전 연구 불가
- 관상동맥은 개구부 스텁만 (관상동맥 유동 해석 불가)
- 심방 중격·난원와, 유스타키오판, 관상정맥동 등 세부 구조 없음
- 방실판 첨판 면적이 실제보다 작음 (건삭이 없어 tenting 단순화)
- 변형은 **규정된 운동학**이지 역학 해석 결과가 아님 (심근 재료모델·전기전도 지연 없음,
  벽 응력 계산 안 됨, 국소 벽운동 이상 표현 불가)
- 대혈관 팽창에 맥파 전파 지연 없음 (모든 위치가 동시에 팽창)

환자 특이적이 필요하면 3D Slicer + vmtk 세그멘테이션, 또는
Zygote Solid 3D Heart / Vascular Model Repository(vascularmodel.com) 사용 권장.

---

## 9. 미해결 / 다음에 할 만한 것

사용자에게 제안했으나 아직 하지 않은 것들:

- [ ] **변형 격자 시계열** — 40 프레임 PLOT3D는 프레임당 41 MB → 전체 1.6 GB.
      현실적으로는 "기준 격자 + 변위장 파일" 형태가 맞음.
      `deform.py`의 변위장을 격자 절점에 그대로 적용 가능하게 설계해 둠
- [ ] **유두근 / 건삭 추가** (LV 와류 정확도에 중요)
- [ ] **좌심방이(LAA) 추가** (심방세동 연구에 필수)
- [ ] 협착·역류 등 **병변 모사** (`params.py`의 해당 혈관 반경 국소 축소)
- [ ] 다른 체격으로 **스케일링** (O_LV, LV_R, LV_C 및 혈관 z값 비율 조정)
- [ ] 솔버별 연결정보 파일 변환 — CFL3D `.inp`, Wind-US `.gbc`, NMF 형식
- [ ] 완전 1-to-1 butterfly 토폴로지 격자 (ICEM/Pointwise 필요)

---

## 10. 시행착오 기록 (같은 실수 반복 방지)

| 문제 | 원인 | 해결 |
|---|---|---|
| 평면 캡에 25 mm 짜리 긴 삼각형 | fan 삼각화 + 아이코스피어 절단 | **Delaunay 캡** + 심실을 loft로 재구성 |
| 심첨 부근 긴 변 | 축방향 균등 샘플링 | **자오선 호길이 균등** 샘플링 |
| 대동맥판 첨판 하나가 유출로로 튐 | 레이캐스팅이 LVOT로 빠져나감 | **해석적 근부 반경 프로파일** 사용 |
| 닫힌 판막이 도메인을 밀봉 못함 | 인접 첨판이 두께 0으로 접촉 | 접합선을 0.8 mm 겹치게 + 반월판 각도 0.07 rad 확장 |
| 대동맥 격자 음의 셀 | 곡률 반경 < 관 반경 (궁 안쪽) | `r ≤ 0.82/κ` 제한 |
| RV 격자 음의 셀 (2~4개) | 심첨에서 초승달 두께 → 0 | θ 범위 3.5% 축소 + 축방향 끝 3.5% 절단 |
| 체적 보정 발산 | 좌/우 solid 오판 (중심점 거리 비교) | `contains(O_LV)` 로 판정 |
| `/mnt/user-data/outputs` 에 numpy `tofile` 실패 | rclone filestore | `fh.write(arr.tobytes())` + 로컬 생성 후 복사 |

---

## 11. 새 세션 시작 프롬프트 예시

> 첨부한 `HANDOFF_심장CFD모델.md` 와 `scripts/` 폴더가 이전 작업 내용입니다.
> `params.py`에 모든 해부학적 치수가, `build.py`의 `FORBIDDEN`에 간섭 검사 쌍이 있습니다.
> 이어서 [ 하고 싶은 일 ] 을 해 주세요.
> 치수를 수정하면 반드시 `python3 build.py lowres` 로
> `clearance violations: none` 을 먼저 확인해 주세요.
