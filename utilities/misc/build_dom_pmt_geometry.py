#!/usr/bin/env python3
"""
Build the real KM3NeT DOM PMT geometry and print absolute positions + directions.

Two independent pieces:
  1. The single-DOM PMT layout (31 PMTs), reproduced from the original survey :
     a handful of MEASURED PMT coordinates define the DOM's
     local axes, and each ring is filled by rotating one measured "seed" PMT in 60 deg
     steps about the cable/heatsink axis (+y). This yields the real, slightly
     aspherical layout (radii ~196-198 mm) with the heatsink at the +y pole.
  2. The DOM centres in the detector (the latest hex A-frame DU geometry).

For every DOM we place the local layout at that centre and print the absolute
(x, y, z) position and unit direction of each PMT. Directions are the DOM-local unit
vectors (all DOMs share the same +y-up orientation), identical at every DOM.

This is the geometry stored in PMTINFO.ratdb as DOMINFO_aframe_spacing5m_hex_latest.
"""

import numpy as np

# ----------------------------------------------------------------------------
# 1. Real single-DOM PMT layout, from the measured survey
# ----------------------------------------------------------------------------

# Measured coordinates (mm) in the raw survey frame. Only a generating subset was
# surveyed: the bottom-pole PMT (#1), one "seed" PMT per ring (#2,8,20,19,31), and
# three more PMTs (#5,11,17) used solely to define the DOM's local axes.
MEAS = {
    1:  (0.0,     0.0,     0.0),      # bottom pole
    2:  (60.25,   85.12,   42.75),    # ring seed
    8:  (185.91,  116.25,  90.14),    # ring seed
    20: (325.54,  35.24,   109.82),   # ring seed
    19: (148.64,  116.98, -11.61),    # ring seed
    31: (304.61,  77.76,   12.40),    # ring seed
    5:  (-4.89,  -107.88, -32.27),    # axis definition only
    11: (34.43,  -77.48,   169.48),   # axis definition only
    17: (128.84,  10.36,  -138.59),   # axis definition only
}
M = {k: np.array(v, float) for k, v in MEAS.items()}
unit = lambda v: v / np.linalg.norm(v)

# DOM-local axes from the hardware: x across the support attachments, z along the
# rope plane, y = z x x  (points from the bottom face up to the cable exit = heatsink).
xhat = unit(M[11] - M[17])
zhat = unit(M[2]  - M[5])
yhat = np.cross(zhat, xhat)

def to_local(k):
    """Measured PMT k expressed in the DOM-local (x, y, z) frame."""
    v = M[k]
    return np.array([v @ xhat, v @ yhat, v @ zhat])

# DOM centre height along y = midpoint of the top and bottom rings.
y_centre = (to_local(31)[1] + to_local(8)[1]) / 2.0

def rot_y(deg):
    """Rotation matrix about the +y (heatsink) axis."""
    t = np.deg2rad(deg); c, s = np.cos(t), np.sin(t)
    return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])

# Fill the 31 PMTs in index order 1..31: the bottom pole, then each ring (one seed
# PMT rotated by 0/60/.../300 deg). (seed_id, rotation_deg) per PMT index:
LAYOUT = {
    1:  (1, 0),
    2:  (2, 0),  3:  (2, 300), 4:  (2, 240), 5:  (2, 180), 6:  (2, 120), 7:  (2, 60),
    8:  (8, 0),  10: (8, 300), 12: (8, 240), 14: (8, 180), 16: (8, 120), 18: (8, 60),
    20: (20, 0), 22: (20, 300),24: (20, 240),26: (20, 180),28: (20, 120),30: (20, 60),
    19: (19, 0), 9:  (19, 300),11: (19, 240),13: (19, 180),15: (19, 120),17: (19, 60),
    31: (31, 0), 21: (31, 300),23: (31, 240),25: (31, 180),27: (31, 120),29: (31, 60),
}

local_pos = np.zeros((31, 3))
for idx in range(1, 32):
    seed, deg = LAYOUT[idx]
    p = to_local(seed) - np.array([0.0, y_centre, 0.0])   # centre on DOM
    local_pos[idx - 1] = rot_y(deg) @ p

# Re-centre on the PMT sphere-centre so the DOM lines up with its glass-sphere centre
# (the real layout is slightly aspherical, so this differs from the ring midpoint by
# ~0.2 mm). Directions are the outward radial unit vectors.
def sphere_centre(P):
    A = np.c_[2 * P, np.ones(len(P))]
    c, *_ = np.linalg.lstsq(A, (P ** 2).sum(1), rcond=None)
    return c[:3]
local_pos -= sphere_centre(local_pos)
local_dir = local_pos / np.linalg.norm(local_pos, axis=1, keepdims=True)

# ----------------------------------------------------------------------------
# 2. DOM centres in the detector (latest hex A-frame DU geometry), mm
#    5 A-frames x 3 layers (y) x 6-DOM hexagon (x-z).
# ----------------------------------------------------------------------------
DOM_CENTERS = [
    (np.float64(1499.9999999999898), np.float64(-2.2282620193436742e-11), np.float64(7198.076211353318)),
    (np.float64(-1500.0000000000136), np.float64(8.185452315956354e-12), np.float64(7198.07621135332)),
    (np.float64(-2999.999999999988), np.float64(-4.160938260611147e-11), np.float64(4600.000000000005)),
    (np.float64(-1499.9999999999932), np.float64(-1.1084466677857563e-10), np.float64(2001.9237886466892)),
    (np.float64(1499.999999999994), np.float64(-4.1313796828035265e-10), np.float64(2001.9237886466904)),
    (np.float64(2999.9999999999936), np.float64(1.6143530956469476e-11), np.float64(4600.000000000004)),
    (np.float64(1499.9999999999952), np.float64(-3000.0000000000055), np.float64(7198.076211353314)),
    (np.float64(-1500.0000000000127), np.float64(-2999.999999999977), np.float64(7198.076211353318)),
    (np.float64(-3000.000000000001), np.float64(-3000.000000000006), np.float64(4599.999999999999)),
    (np.float64(-1499.999999999992), np.float64(-3000.0000000000036), np.float64(2001.9237886466833)),
    (np.float64(1500.0000000000098), np.float64(-3000.0000000000064), np.float64(2001.9237886466767)),
    (np.float64(3000.000000000003), np.float64(-2999.999999999999), np.float64(4600.000000000004)),
    (np.float64(1499.9999999999875), np.float64(-6000.00000000001), np.float64(7198.07621135331)),
    (np.float64(-1499.9999999999764), np.float64(-6000.000000000015), np.float64(7198.076211353303)),
    (np.float64(-3000.0000000000005), np.float64(-5999.999999999996), np.float64(4600.000000000005)),
    (np.float64(-1499.9999999999916), np.float64(-6000.000000000002), np.float64(2001.9237886466806)),
    (np.float64(1500.0000000000234), np.float64(-6000.000000000006), np.float64(2001.9237886466813)),
    (np.float64(3000.0000000000055), np.float64(-5999.999999999982), np.float64(4600.000000000021)),
    (np.float64(1499.9999999428876), np.float64(2.372975359321572e-06), np.float64(18198.07621135679)),
    (np.float64(-1500.0000000001737), np.float64(1.4551915228366852e-11), np.float64(18198.076211353306)),
    (np.float64(-3000.000000000044), np.float64(2.7284841053187847e-11), np.float64(15600.000000000033)),
    (np.float64(-1499.9999999999177), np.float64(-4.001776687800884e-11), np.float64(13001.923788646669)),
    (np.float64(1499.9999999998317), np.float64(1.9099388737231493e-11), np.float64(13001.92378864672)),
    (np.float64(3000.0000000002115), np.float64(-1.3733369996771216e-10), np.float64(15599.999999999929)),
    (np.float64(1499.9999998669196), np.float64(-2999.9999968321463), np.float64(18198.076211884836)),
    (np.float64(-1500.0000000000825), np.float64(-2999.999999999978), np.float64(18198.076211353364)),
    (np.float64(-3000.0000000000305), np.float64(-3000.000000000018), np.float64(15600.000000000022)),
    (np.float64(-1499.9999999999918), np.float64(-3000.000000000018), np.float64(13001.923788646658)),
    (np.float64(1499.9999999998875), np.float64(-2999.9999999999372), np.float64(13001.923788646714)),
    (np.float64(3000.0000000002474), np.float64(-3000.000000000089), np.float64(15599.999999999949)),
    (np.float64(1499.9999998342616), np.float64(-5999.999997424742), np.float64(18198.0762122147)),
    (np.float64(-1499.9999999999645), np.float64(-6000.000000000024), np.float64(18198.076211353306)),
    (np.float64(-2999.999999999998), np.float64(-6000.000000000003), np.float64(15600.000000000013)),
    (np.float64(-1499.9999999999097), np.float64(-6000.000000000013), np.float64(13001.92378864664)),
    (np.float64(1499.9999999999002), np.float64(-5999.999999999969), np.float64(13001.923788646694)),
    (np.float64(3000.000000000209), np.float64(-6000.000000000091), np.float64(15599.999999999973)),
    (np.float64(1499.9999999998108), np.float64(5.093170329928398e-11), np.float64(29198.076211353276)),
    (np.float64(-1500.000000000059), np.float64(-7.09405867382884e-11), np.float64(29198.076211353145)),
    (np.float64(-3000.0000000001114), np.float64(5.275069270282984e-11), np.float64(26600.00000000003)),
    (np.float64(-1499.9999999998959), np.float64(1.8189894035458565e-11), np.float64(24001.923788646658)),
    (np.float64(1499.999999999912), np.float64(-8.913048077374697e-11), np.float64(24001.923788646676)),
    (np.float64(3000.000000000224), np.float64(1.8189894035458565e-11), np.float64(26599.999999999847)),
    (np.float64(1500.0000000001337), np.float64(-3000.00000000004), np.float64(29198.076211353284)),
    (np.float64(-1499.9999999998938), np.float64(-3000.0000000000673), np.float64(29198.076211353262)),
    (np.float64(-3000.0000000002165), np.float64(-2999.9999999999272), np.float64(26600.000000000076)),
    (np.float64(-1500.0000000000314), np.float64(-2999.999999999922), np.float64(24001.92378864668)),
    (np.float64(1499.9999999994805), np.float64(-2999.9999999999654), np.float64(24001.923788646734)),
    (np.float64(3000.000000000351), np.float64(-3000.0000000001273), np.float64(26599.99999999983)),
    (np.float64(1500.0000003828775), np.float64(-6000.000006843966), np.float64(29198.076209929626)),
    (np.float64(-1499.9999999998872), np.float64(-6000.0), np.float64(29198.076211353247)),
    (np.float64(-2999.9999999999845), np.float64(-6000.000000000047), np.float64(26599.999999999978)),
    (np.float64(-1499.9999999999445), np.float64(-5999.999999999944), np.float64(24001.923788646676)),
    (np.float64(1500.000000228351), np.float64(-6000.000003845391), np.float64(24001.92378767279)),
    (np.float64(3000.00000000012), np.float64(-5999.999999999947), np.float64(26599.99999999986)),
    (np.float64(1499.9999999994925), np.float64(6.293703336268663e-10), np.float64(40198.07621135379)),
    (np.float64(-1499.9999999992322), np.float64(2.9467628337442875e-10), np.float64(40198.07621135319)),
    (np.float64(-3000.000000000412), np.float64(1.4551915228366852e-10), np.float64(37600.00000000029)),
    (np.float64(-1500.0000000004875), np.float64(5.384208634495735e-10), np.float64(35001.92378864708)),
    (np.float64(1499.9999999997672), np.float64(-1.5643308870494366e-10), np.float64(35001.923788646745)),
    (np.float64(3000.000000012499), np.float64(-4.7607318265363574e-07), np.float64(37599.99999999875)),
    (np.float64(1500.000000534592), np.float64(-3000.0000137370444), np.float64(40198.076210311745)),
    (np.float64(-1499.999999999906), np.float64(-2999.999999999731), np.float64(40198.07621135332)),
    (np.float64(-3000.000000000428), np.float64(-2999.99999999952), np.float64(37600.00000000019)),
    (np.float64(-1500.0000000004297), np.float64(-2999.9999999997963), np.float64(35001.92378864708)),
    (np.float64(1500.0000005255397), np.float64(-3000.0000132398454), np.float64(35001.9237874931)),
    (np.float64(2999.9999999907945), np.float64(-2999.9999996825936), np.float64(37600.00000002548)),
    (np.float64(1500.000000479037), np.float64(-6000.0000094389725), np.float64(40198.076209929044)),
    (np.float64(-1500.0000000000357), np.float64(-5999.999999999767), np.float64(40198.07621135336)),
    (np.float64(-3000.000000000639), np.float64(-5999.999999999865), np.float64(37600.00000000023)),
    (np.float64(-1500.0000000007815), np.float64(-5999.9999999998), np.float64(35001.92378864691)),
    (np.float64(1500.000000456073), np.float64(-6000.000008620926), np.float64(35001.92378715178)),
    (np.float64(2999.999999972052), np.float64(-5999.999999305339), np.float64(37600.000000112675)),
    (np.float64(1500.0000005347267), np.float64(-1.8965874915011227e-05), np.float64(51198.076211341424)),
    (np.float64(-1500.0000000001853), np.float64(-1.0913936421275139e-10), np.float64(51198.07621135334)),
    (np.float64(-3000.0000000015725), np.float64(1.1568772606551647e-09), np.float64(48600.00000000093)),
    (np.float64(-1499.9999999996405), np.float64(5.638867150992155e-10), np.float64(46001.92378864692)),
    (np.float64(1500.0000005590844), np.float64(-1.995550701394677e-05), np.float64(46001.923788632805)),
    (np.float64(3000.000000211205), np.float64(-8.179609721992165e-06), np.float64(48599.999999988475)),
    (np.float64(1500.0000005075237), np.float64(-3000.0000135184346), np.float64(51198.076210549385)),
    (np.float64(-1499.9999989291264), np.float64(-3000.0000235941843), np.float64(51198.0762100067)),
    (np.float64(-3000.0000000007535), np.float64(-2999.9999999996544), np.float64(48600.000000000386)),
    (np.float64(-1500.000000001302), np.float64(-2999.9999999994434), np.float64(46001.923788646935)),
    (np.float64(1500.00000052662), np.float64(-3000.0000137516545), np.float64(46001.92378773545)),
    (np.float64(3000.000000190955), np.float64(-3000.0000056917816), np.float64(48599.99999963829)),
    (np.float64(1500.0000004697636), np.float64(-6000.000009999709), np.float64(51198.07621017028)),
    (np.float64(-1499.9999989924183), np.float64(-6000.000017278708), np.float64(51198.07620936121)),
    (np.float64(-2999.9999986581097), np.float64(-6000.000020485997), np.float64(48599.99999755844)),
    (np.float64(-1499.9999988654513), np.float64(-6000.00001833764), np.float64(46001.92378629591)),
    (np.float64(1500.000000481081), np.float64(-6000.000009870266), np.float64(46001.9237873457)),
    (np.float64(3000.000000167989), np.float64(-6000.000004037858), np.float64(48599.99999949195)),
]

# ----------------------------------------------------------------------------
# 3. Place the DOM at every centre; print absolute positions and directions
# ----------------------------------------------------------------------------
print(f"# {len(DOM_CENTERS)} DOMs x 31 PMTs = {len(DOM_CENTERS)*31} PMTs")
print(f"# {'idx':>5} {'dom':>4} {'pmt':>3}  "
      f"{'x':>11} {'y':>11} {'z':>11}   {'dir_x':>9} {'dir_y':>9} {'dir_z':>9}")
idx = 0
for dom, centre in enumerate(DOM_CENTERS):
    c = np.array(centre, float)
    for j in range(31):
        p = c + local_pos[j]
        d = local_dir[j]
        print(f"  {idx:>5} {dom:>4} {j:>3}  "
              f"{p[0]:>11.3f} {p[1]:>11.3f} {p[2]:>11.3f}   "
              f"{d[0]:>9.5f} {d[1]:>9.5f} {d[2]:>9.5f}")
        idx += 1
