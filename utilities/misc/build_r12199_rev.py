#!/usr/bin/env python3
"""
Generate (and plot) the profile for the r12199_rev PMT model — RAT "revolution" construction.

Built faithfully from the Hamamatsu R12199 datasheet (TPMH1356E, Fig. 4):

    Ø80 mm    hemispherical bulb (outer glass)
    Ø72 mm    minimum photocathode effective area
    Ø51.9 mm  cylindrical neck / tube
    Ø56.5 mm  base
    110 mm    overall body length (photocathode apex -> base bottom)
    ~3 mm     glass wall thickness

A "revolution" PMT is two polylines in the (rho, z) plane, spun around the z-axis:
the OUTER glass envelope and the INNER vacuum surface. `build_profiles()` returns both;
`emit_ratdb()` prints them for PMT.ratdb; `plot()` draws the two-panel figure.

Coordinate / placement convention
----------------------------------
z is the PMT axis; +z points OUT toward the light (toward the photocathode).
We put the dome APEX at z = +4 mm, matching the old toroidal r12199, so that when a PMT
is placed at the existing DOM mount positions (~196-198 mm from the DOM centre) the
photocathode still lands at the inner glass (~203 mm). The bulb is a hemisphere of
radius 40 mm, so its centre / widest ring (the "equator") sits 40 mm below the apex at
z = -36 mm.

NOTE: dynode dimensions are NOT computed here and are NOT from the datasheet (which gives
no dynode mechanics). They are placeholders; the dynode is an opaque internal structure,
optically irrelevant for a Cherenkov simulation.
"""

import math

# --------------------------------------------------------------------------
# 1. Datasheet dimensions (mm)
# --------------------------------------------------------------------------
BULB_RADIUS = 40.0    # Ø80 bulb
NECK_RADIUS = 25.95   # Ø51.9 tube
BASE_RADIUS = 28.25   # Ø56.5 base
WALL        = 3.0     # glass wall thickness
BODY_LENGTH = 110.0   # apex -> base bottom

# --------------------------------------------------------------------------
# 2. Local-frame landmarks along the axis (mm), from the convention above
# --------------------------------------------------------------------------
Z_APEX     = 4.0                    # dome tip (front, toward the light)
Z_EQUATOR  = Z_APEX - BULB_RADIUS   # -36 : hemisphere centre = widest ring
Z_NECK     = -51.0                  # shoulder has narrowed to the tube by here
Z_BASE     = Z_APEX - BODY_LENGTH   # -106 : base bottom (110 mm below the apex)
Z_TUBE_END = Z_BASE + 4.0           # -102 : tube ends, then a short flare to the base


# --------------------------------------------------------------------------
# 3. Profile-building helpers
# --------------------------------------------------------------------------
def dome_arc(radius, n=40):
    """Hemisphere from apex to equator as (rho, z) points.
    Centred on the axis at z = Z_EQUATOR; phi is the polar angle from the apex."""
    pts = []
    for k in range(n):
        phi = math.radians(2.0 + (90.0 - 2.0) * k / (n - 1))
        pts.append((radius * math.sin(phi), Z_EQUATOR + radius * math.cos(phi)))
    return pts


def straight(p0, p1, n=4):
    """Straight run from p0 to p1 (excluding p0) as (rho, z) points."""
    return [(p0[0] + (p1[0] - p0[0]) * i / (n - 1),
             p0[1] + (p1[1] - p0[1]) * i / (n - 1)) for i in range(1, n)]


def smooth_shoulder(rho_top, rho_bot, z_top, z_bot, n=16):
    """Smoothly narrow the glass from the bulb equator (rho_top) to the tube (rho_bot),
    over z_top->z_bot. Uses a smoothstep (3t^2-2t^3) so the curve is tangent to the
    hemisphere at the top and to the cylinder at the bottom — no sharp kinks, matching
    the datasheet's continuously-curved envelope."""
    pts = []
    for i in range(1, n):
        t = i / (n - 1)
        s = 3 * t * t - 2 * t * t * t          # 0 -> 1, zero slope at both ends
        pts.append((rho_top + (rho_bot - rho_top) * s, z_top + (z_bot - z_top) * t))
    return pts


def build_profiles():
    """Return (outer, inner) profiles, each a list of (rho, z) points, front to back."""
    # OUTER glass:  dome -> shoulder -> tube -> base
    outer  = dome_arc(BULB_RADIUS)
    outer += smooth_shoulder(BULB_RADIUS, NECK_RADIUS, Z_EQUATOR, Z_NECK)   # curved shoulder Ø80->Ø51.9
    outer += straight((NECK_RADIUS, Z_NECK),    (NECK_RADIUS, Z_TUBE_END))  # cylindrical tube
    outer += [(BASE_RADIUS, Z_TUBE_END - 1.0), (BASE_RADIUS, Z_BASE)]       # short flare to Ø56.5 base

    # INNER vacuum: same profile shrunk inward by the wall thickness
    inner  = dome_arc(BULB_RADIUS - WALL)
    inner += smooth_shoulder(BULB_RADIUS - WALL, NECK_RADIUS - WALL, Z_EQUATOR, Z_NECK)
    inner += straight((NECK_RADIUS - WALL, Z_NECK),    (NECK_RADIUS - WALL, Z_TUBE_END))
    inner += [(NECK_RADIUS - WALL, Z_TUBE_END - 1.0)]
    return outer, inner


# --------------------------------------------------------------------------
# 4. Emit arrays for PMT.ratdb
# --------------------------------------------------------------------------
def emit_ratdb(outer, inner):
    def line(name, vals):
        return '"%s": [ %s ],' % (name, ", ".join("%.3f" % v for v in vals))
    print('// --- paste under index "r12199_rev" in PMT.ratdb ---')
    print(line("rho_inner", [p[0] for p in inner]))
    print(line("z_inner",   [p[1] for p in inner]))
    print(line("rho_edge",  [p[0] for p in outer]))
    print(line("z_edge",    [p[1] for p in outer]))
    print()
    print("// dynode fields are placeholders (not datasheet, optically irrelevant):")
    print('"dynode_radius": 20.0,\n"dynode_top": -46.0,\n"dynode_height": 50.0,')


# --------------------------------------------------------------------------
# 5. Plot: (A) the profile, (B) the PMT seated in the DOM glass
# --------------------------------------------------------------------------
def plot(outer, inner, mount_radius=199.0, outfile="r12199_rev_geometry.png"):
    """mount_radius = radial distance from the DOM centre to the PMT's local origin;
    with the apex at z=+4 this puts the photocathode at mount_radius+4 (~the inner glass)."""
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    oz  = [p[1] for p in outer]; orho = [p[0] for p in outer]
    iz  = [p[1] for p in inner]; irho = [p[0] for p in inner]
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(15, 7))

    # Panel A: the profile (x = z, y = rho)
    axA.plot(oz, orho, '-', color='#1e8449', lw=2)
    axA.plot(oz, [-r for r in orho], '-', color='#1e8449', lw=2)
    axA.fill_between(oz, [-r for r in orho], orho, color='#27ae60', alpha=0.10)
    axA.plot(iz, irho, '--', color='#7f8c8d', lw=1, label='inner (vacuum) surface')
    axA.plot(iz, [-r for r in irho], '--', color='#7f8c8d', lw=1)
    axA.annotate("Ø80 curved dome", (-15, 37), xytext=(-15, 58), fontsize=10,
                 color='#1e8449', ha='center', arrowprops=dict(arrowstyle='->', color='#1e8449'))
    axA.annotate("Ø51.9 tube", (-80, NECK_RADIUS), xytext=(-80, 42), fontsize=9,
                 ha='center', color='#555', arrowprops=dict(arrowstyle='->', color='#555'))
    axA.annotate("Ø56.5 base", (Z_BASE, BASE_RADIUS), xytext=(-95, 12), fontsize=9,
                 ha='center', color='#555', arrowprops=dict(arrowstyle='->', color='#555'))
    axA.axhline(0, color='k', lw=0.6, ls=':')
    axA.set_title("r12199_rev — Ø80 dome, Ø51.9 tube, 110 mm (datasheet Fig. 4)", fontsize=11)
    axA.set_xlabel("z  (mm, PMT axis — front = +z)"); axA.set_ylabel("ρ  (mm)")
    axA.set_aspect('equal'); axA.grid(alpha=0.25); axA.legend(loc='lower left', fontsize=8)

    # Panel B: seated in the DOM (radial distance from DOM centre = mount_radius + z)
    ox = [mount_radius + z for z in oz]
    axB.plot(ox, orho, '-', color='#1e8449', lw=2, label='r12199_rev (real dome)')
    axB.plot(ox, [-r for r in orho], '-', color='#1e8449', lw=2)
    axB.fill_between(ox, [-r for r in orho], orho, color='#27ae60', alpha=0.12)
    th = np.linspace(-0.45, 0.45, 200)
    for R, ls, lab in [(203, '-', 'inner glass (203 mm)'), (216, '--', 'outer glass (216 mm)')]:
        axB.plot(R * np.cos(th), R * np.sin(th), ls, color='#2e86c1', lw=1.6, label=lab)
    axB.set_title("Seated at mount %.0f mm: dome nestles inside the inner glass" % mount_radius, fontsize=11)
    axB.set_xlabel("radial distance from DOM centre (mm)"); axB.set_ylabel("lateral (mm)")
    axB.set_aspect('equal'); axB.grid(alpha=0.25); axB.legend(loc='lower left', fontsize=8)
    axB.set_xlim(80, 235); axB.set_ylim(-70, 70)

    plt.tight_layout(); plt.savefig(outfile, dpi=130, bbox_inches='tight'); plt.close()
    return outfile


# --------------------------------------------------------------------------
if __name__ == "__main__":
    outer, inner = build_profiles()
    emit_ratdb(outer, inner)
    out = plot(outer, inner)
    print("\n// summary: bulb Ø%.1f, tube Ø%.1f, base Ø%.1f, length %.1f mm  ->  wrote %s"
          % (2 * max(p[0] for p in outer), 2 * NECK_RADIUS, 2 * BASE_RADIUS,
             max(p[1] for p in outer) - min(p[1] for p in outer), out))
