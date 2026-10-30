"""
params.py -- anatomical dimensions for a 'standard' adult (male, 175 cm, ~70 kg, BSA 1.85 m2)
All lengths in mm.  Body frame:  +x = LEFT, +y = ANTERIOR, +z = SUPERIOR (height above floor)
"""
import numpy as np

# ======================================================================
# GLOBAL CARDIAC AXIS
# ======================================================================
D_LV = np.array([0.501, 0.301, -0.812])          # base -> apex direction
D_LV = D_LV / np.linalg.norm(D_LV)

O_LV = np.array([46.0, 14.0, 1290.0])            # LV cavity ellipsoid centre
LV_R = 26.0                                       # LV endocardial short semi-axis
LV_C = 52.0                                       # LV endocardial long semi-axis
LV_BASE_OFF = 32.0                                # truncation distance from centre toward base
LV_WALL = 10.0                                    # LV wall thickness (mid-wall, diastole)

B_LV = O_LV - LV_BASE_OFF * D_LV                  # LV base-plane centre
A_LV = O_LV + LV_C * D_LV                         # LV endocardial apex

# unit vector perpendicular to the LV axis, pointing anterior-right  (RV side)
_w = np.array([-0.75, 0.66, 0.0])
_w = _w - np.dot(_w, D_LV) * D_LV
W_RV = _w / np.linalg.norm(_w)

# RV outer envelope
RV_SHIFT = 20.0
O_RV = O_LV + RV_SHIFT * W_RV
RV_R = 41.0
RV_C = 54.0
RV_BASE_OFF = 30.0
B_RV = O_LV - RV_BASE_OFF * D_LV

# ======================================================================
# ATRIA  (ellipsoids: centre, semi-axes)
# ======================================================================
LA_C = np.array([26.0, -26.0, 1350.0]);  LA_R = np.array([25.0, 17.0, 26.0])
RA_C = np.array([-33.0,   6.0, 1330.0]); RA_R = np.array([21.0, 22.0, 27.0])

# ======================================================================
# LEFT-HEART OUTFLOW : LVOT -> aortic root -> ascending -> arch -> descending
#  (first points lie INSIDE the LV cavity so the union is clean)
# ======================================================================
AORTA_PTS = [
    (37.0,  12.0, 1305.0),   # intracavitary LVOT
    (22.0,   9.0, 1318.0),   # LVOT
    (12.0,   7.0, 1326.0),   # subvalvular
    ( 5.0,   6.0, 1332.0),   # AORTIC ANNULUS
    ( 3.0,   6.0, 1342.0),   # sinuses of Valsalva (max)
    (-1.0,   5.0, 1354.0),   # sino-tubular junction
    (-12.0,  3.0, 1370.0),   # ascending
    (-19.0,  2.0, 1392.0),   # upper ascending
    (-18.0, -6.0, 1410.0),   # proximal arch
    (-10.0,-22.0, 1421.0),   # arch apex
    (  0.0,-36.0, 1417.0),   # distal arch
    (  8.0,-48.0, 1400.0),   # isthmus
    ( 11.0,-55.0, 1372.0),   # proximal descending
    ( 12.0,-58.0, 1320.0),
    ( 12.0,-58.0, 1250.0),
    ( 12.0,-58.0, 1180.0),   # OUTLET (diaphragmatic level)
]
AORTA_R = [21.0, 15.5, 13.2, 12.5, 16.8, 14.2, 15.0, 14.6,
           14.2, 13.6, 12.8, 11.8, 11.5, 11.2, 10.9, 10.6]
AO_ANNULUS_IDX = 3
AO_SINUS_IDX = 4
AO_STJ_IDX = 5

# aortic arch branches: (control points, radii)
BRANCHES = {
    "brachiocephalic_trunk": ([(-17.0, -3.0, 1405.0), (-20.0, -8.0, 1428.0),
                               (-26.0,-12.0, 1452.0), (-31.0,-14.0, 1474.0)],
                              [7.5, 6.8, 6.3, 6.0]),
    "left_common_carotid":   ([(-10.0,-21.0, 1418.0), (-9.0, -25.0, 1440.0),
                               (-8.0, -27.0, 1462.0), (-7.0, -28.0, 1480.0)],
                              [4.8, 4.3, 4.1, 4.0]),
    "left_subclavian":       ([(-1.0, -34.0, 1416.0), (4.0, -38.0, 1438.0),
                               (12.0, -41.0, 1458.0), (20.0, -43.0, 1474.0)],
                              [5.4, 5.0, 4.7, 4.5]),
}

# coronary ostia stubs (short outlets off the sinuses)
CORONARIES = {
    "left_coronary_ostium":  ([(3.0, 6.0, 1343.0), (13.0, 12.0, 1341.0), (21.0, 15.0, 1338.0)],
                              [2.6, 2.3, 2.2]),
    "right_coronary_ostium": ([(3.0, 6.0, 1341.0), (-4.0, 16.0, 1343.0), (-9.0, 23.0, 1344.0)],
                              [2.3, 2.0, 1.9]),
}

# ======================================================================
# MITRAL INFLOW TRACT  (LA -> LV)
# ======================================================================
MITRAL_PTS = [(26.0, -20.0, 1344.0), (28.0, -8.0, 1334.0),
              (30.0,  2.0, 1322.0),  (33.0, 12.0, 1306.0)]
MITRAL_R   = [21.0, 17.0, 16.0, 20.0]
MV_ANNULUS_IDX = 2            # annulus plane
MV_ANNULUS_R = 16.0

# pulmonary veins : (control points, radii)   -- inlets of the left circuit
PULM_VEINS = {
    "pv_right_superior": ([(16.0,-30.0,1362.0), (2.0,-35.0,1363.0), (-14.0,-40.0,1364.0), (-32.0,-43.0,1364.0)],
                          [8.0, 6.8, 6.4, 6.2]),
    "pv_right_inferior": ([(16.0,-30.0,1336.0), (2.0,-35.0,1334.0), (-14.0,-40.0,1332.0), (-32.0,-43.0,1330.0)],
                          [7.6, 6.4, 6.0, 5.8]),
    "pv_left_superior":  ([(36.0,-30.0,1362.0), (48.0,-42.0,1372.0), (60.0,-50.0,1378.0), (72.0,-54.0,1382.0)],
                          [7.6, 6.4, 6.1, 5.9]),
    "pv_left_inferior":  ([(36.0,-30.0,1336.0), (48.0,-42.0,1330.0), (60.0,-50.0,1326.0), (72.0,-54.0,1322.0)],
                          [7.2, 6.0, 5.7, 5.5]),
}

# ======================================================================
# RIGHT-HEART OUTFLOW : RVOT -> pulmonary valve -> MPA -> LPA / RPA
# ======================================================================
RVOT_START = O_LV - 20.0 * D_LV + 44.0 * W_RV      # inside the RV crescent
PA_PTS = [tuple(RVOT_START),
          (11.0, 39.0, 1322.0),   # infundibulum
          (13.0, 35.0, 1341.0),
          (18.0, 28.0, 1356.0),   # PULMONARY ANNULUS
          (19.0, 24.0, 1366.0),   # pulmonary sinuses
          (18.0, 14.0, 1377.0),   # main pulmonary artery
          (14.0,  2.0, 1385.0),
          (11.0,-11.0, 1390.0)]   # bifurcation
PA_R = [19.0, 14.0, 12.2, 11.5, 13.2, 13.6, 13.4, 13.0]
PV_ANNULUS_IDX = 3
PA_BIF = np.array([11.0, -11.0, 1390.0])

RPA_PTS = [(11.0,-11.0,1390.0), (-6.0,-26.0,1388.0), (-26.0,-33.0,1385.0),
           (-44.0,-35.0,1381.0), (-60.0,-35.0,1377.0)]
RPA_R   = [12.5, 11.0, 10.8, 10.6, 10.4]

LPA_PTS = [(11.0,-11.0,1390.0), (26.0,-22.0,1392.0), (42.0,-28.0,1392.0),
           (58.0,-30.0,1390.0)]
LPA_R   = [12.0, 10.2, 9.9, 9.7]

# ======================================================================
# TRICUSPID INFLOW  (RA -> RV)   and caval veins (inlets of right circuit)
# ======================================================================
TRICUSPID_PTS = [(-31.0, 10.0, 1326.0), (-24.0, 14.0, 1320.0),
                 (-13.0, 24.0, 1308.0), (6.0, 35.0, 1288.0)]
TRICUSPID_R   = [24.0, 18.5, 17.5, 19.0]
TV_ANNULUS_IDX = 2
TV_ANNULUS_R  = 18.0

CAVAE = {
    "svc": ([(-40.0, 8.0, 1345.0), (-43.0, 6.0, 1385.0), (-42.0, 3.0, 1420.0), (-40.0, 1.0, 1452.0)],
            [12.0, 9.8, 9.6, 9.5]),
    "ivc": ([(-33.0, 0.0, 1312.0), (-31.0, -6.0, 1285.0), (-28.0,-14.0, 1255.0), (-26.0,-20.0, 1228.0)],
            [14.0, 12.2, 11.8, 11.5]),
}

# ======================================================================
# MESH RESOLUTION
# ======================================================================
TARGET_EDGE = 0.75        # mm  -- target triangle edge length
