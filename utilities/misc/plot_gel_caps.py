#!/usr/bin/env python3
"""
Explanatory plot for the per-PMT gel-cap coupling model (NOT a generator).

Shows, in the DOM cross-section, how a gel cap built as a RAT `sphere`
(r_min / r_max / theta_delta, rotated onto the PMT axis) seats between the PMT
photocathode dome and the inner glass — the primitive-only alternative to a
Boolean meniscus.

The PMT profile is imported from build_r12199_rev.py so the shape drawn here is
exactly the simulated r12199_rev. Iterate on GEL_THICK / MOUNT_RECESSED /
PHOTOCATHODE_HALF_ANGLE and re-run to see the fit against the glass and between
neighbouring caps before committing to geo entries.

    python3 plot_gel_caps.py   ->  gel_caps_geometry.png
"""

import os
import sys
import math

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_r12199_rev import build_profiles, BULB_RADIUS, Z_EQUATOR

# --------------------------------------------------------------------------
# DOM + gel parameters (mm, deg) — the knobs to iterate
# --------------------------------------------------------------------------
INNER_GLASS = 203.0
OUTER_GLASS = 216.0
GEL_THICK   = 3.0               # r_max - r_min of the cap shell
MOUNT_RECESSED = 196.0          # 3 mm inside the touching-dome mount (199) so the
                                # cap apex lands on the glass, not past it
PHOTOCATHODE_HALF_ANGLE = 64.0  # asin(36/40): Ø72 photocathode on the Ø80 dome

R_IN  = BULB_RADIUS             # 40 : cap inner face hugs the dome outer surface
R_OUT = BULB_RADIUS + GEL_THICK # 43 : + gel thickness

GREEN, BLUE, GOLD, AIR = '#1e8449', '#2e86c1', '#d4ac0d', '#aed6f1'


# --------------------------------------------------------------------------
# Geometry helpers — everything in the DOM cross-section (x = radial, y = lateral)
# --------------------------------------------------------------------------
def frame(phi_deg):
    """Unit axis (outward) and lateral vectors for a PMT pointing at phi."""
    p = math.radians(phi_deg)
    return (np.array([math.cos(p), math.sin(p)]),
            np.array([-math.sin(p), math.cos(p)]))


def place_profile(profile, mount_radius, phi_deg):
    """Local (rho, z) profile -> the two world edges of the PMT body."""
    ax, perp = frame(phi_deg)
    up   = np.array([(mount_radius + z) * ax + rho * perp for rho, z in profile])
    down = np.array([(mount_radius + z) * ax - rho * perp for rho, z in profile])
    return up, down


def dome_centre(mount_radius, phi_deg):
    ax, _ = frame(phi_deg)
    return (mount_radius + Z_EQUATOR) * ax          # hemisphere centre = local z=-36


def cap_shell(C, phi_deg):
    """Closed polygon of the gel cap: inner arc (on the dome) + outer arc, both
    swept over +/-PHOTOCATHODE_HALF_ANGLE about the PMT axis."""
    ax, perp = frame(phi_deg)
    a = np.linspace(-math.radians(PHOTOCATHODE_HALF_ANGLE),
                    math.radians(PHOTOCATHODE_HALF_ANGLE), 80)
    inner = np.array([C + R_IN  * (math.cos(t) * ax + math.sin(t) * perp) for t in a])
    outer = np.array([C + R_OUT * (math.cos(t) * ax + math.sin(t) * perp) for t in a])
    return np.vstack([inner, outer[::-1]]), inner, outer


def glass_arc(R, phi_deg, span=42):
    t = np.radians(np.linspace(phi_deg - span, phi_deg + span, 300))
    return R * np.cos(t), R * np.sin(t)


def fill_interior(ax, x_left=110):
    """Light wash for the air-filled dom_inner (everything inside the glass)."""
    gx, gy = glass_arc(INNER_GLASS, 0)
    ax.fill(np.r_[gx, x_left, x_left], np.r_[gy, gy[-1], gy[0]],
            color=AIR, alpha=0.18, zorder=0)


# --------------------------------------------------------------------------
# Plot
# --------------------------------------------------------------------------
def main(outfile="gel_caps_geometry.png"):
    outer_prof, _ = build_profiles()
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(15, 7))

    # ---- Panel A: one PMT, the cap, and the air it leaves ----
    fill_interior(axA)
    axA.plot([], [], color=AIR, lw=6, alpha=0.6, label='air (dom_inner)')
    for R, ls, lab in [(INNER_GLASS, '-', 'inner glass (203)'),
                       (OUTER_GLASS, '--', 'outer glass (216)')]:
        axA.plot(*glass_arc(R, 0), ls, color=BLUE, lw=1.8, label=lab)

    up, down = place_profile(outer_prof, MOUNT_RECESSED, 0)
    axA.plot(up[:, 0], up[:, 1], '-', color=GREEN, lw=2, label='r12199_rev PMT')
    axA.plot(down[:, 0], down[:, 1], '-', color=GREEN, lw=2)
    axA.fill(np.concatenate([up[:, 0], down[::-1, 0]]),
             np.concatenate([up[:, 1], down[::-1, 1]]), color=GREEN, alpha=0.10)

    C = dome_centre(MOUNT_RECESSED, 0)
    poly, inner, outer = cap_shell(C, 0)
    axA.fill(poly[:, 0], poly[:, 1], color=GOLD, alpha=0.6, zorder=5,
             label='gel cap (borosilicate, n≈1.47)')

    axA.plot(*C, 'k+', ms=9); axA.annotate('dome centre C', C, xytext=(120, -18),
             fontsize=8, arrowprops=dict(arrowstyle='->'))
    axA.annotate('cap apex touches\nglass (203)', (203, 0), xytext=(207, 20),
                 fontsize=8, color=GOLD, ha='left',
                 arrowprops=dict(arrowstyle='->', color=GOLD))
    axA.annotate('r_min=40 hugs dome\nr_max=43 (+3 mm gel)\nθ ≤ 64° (Ø72 cathode)',
                 (172, 22), xytext=(120, 44), fontsize=8, color='#7d6608',
                 arrowprops=dict(arrowstyle='->', color=GOLD))
    axA.annotate('AIR gap at the cathode edge\n(ideal gel / Boolean\nwould fill this to glass)',
                 (188, -37), xytext=(200, -60), fontsize=8, color='#1b4f72', ha='left',
                 arrowprops=dict(arrowstyle='->', color=BLUE))

    axA.set_title('Gel cap = spherical shell coating the photocathode', fontsize=11)
    axA.set_xlabel('radial distance from DOM centre (mm)'); axA.set_ylabel('lateral (mm)')
    axA.set_aspect('equal'); axA.grid(alpha=0.25); axA.legend(loc='lower left', fontsize=8)
    axA.set_xlim(110, 232); axA.set_ylim(-70, 70)

    # ---- Panel B: two neighbours — caps are siblings, must not overlap ----
    fill_interior(axB)
    axB.plot(*glass_arc(INNER_GLASS, 0), '-', color=BLUE, lw=1.8, label='inner glass (203)')
    for phi in (+14, -14):
        up, down = place_profile(outer_prof, MOUNT_RECESSED, phi)
        axB.plot(up[:, 0], up[:, 1], '-', color=GREEN, lw=1.8)
        axB.plot(down[:, 0], down[:, 1], '-', color=GREEN, lw=1.8)
        axB.fill(np.concatenate([up[:, 0], down[::-1, 0]]),
                 np.concatenate([up[:, 1], down[::-1, 1]]), color=GREEN, alpha=0.10)
        C = dome_centre(MOUNT_RECESSED, phi)
        poly, _, _ = cap_shell(C, phi)
        axB.fill(poly[:, 0], poly[:, 1], color=GOLD, alpha=0.55, zorder=5)
    axB.plot([], [], color=GOLD, lw=6, alpha=0.55, label='gel caps (siblings)')

    axB.annotate('adjacent caps must clear\neach other AND the PMTs\n(overlap checker enforces)',
                 (196, 0), xytext=(150, 40), fontsize=8, ha='center',
                 arrowprops=dict(arrowstyle='->'))
    axB.set_title('Two neighbours: caps as non-overlapping siblings', fontsize=11)
    axB.set_xlabel('radial distance from DOM centre (mm)'); axB.set_ylabel('lateral (mm)')
    axB.set_aspect('equal'); axB.grid(alpha=0.25); axB.legend(loc='lower left', fontsize=8)
    axB.set_xlim(110, 232); axB.set_ylim(-70, 70)

    plt.tight_layout(); plt.savefig(outfile, dpi=130, bbox_inches='tight'); plt.close()
    return outfile


def plot_fat_shell(outfile="fat_shell_geometry.png"):
    """The 'fat meniscus' interpretation: gel fills the band r=R_SHELL_IN..203 up
    to the dome faces (= shell minus domes). Drawn by laying down the full band and
    overdrawing the two domes to carve them out — which is exactly the Boolean."""
    outer_prof, _ = build_profiles()
    fig, ax = plt.subplots(figsize=(9.5, 8))

    MOUNT = 199.0          # dome apex at 203: touches the glass (not recessed)
    R_SHELL_IN = 187.0     # inner bound of the shell = photocathode-edge depth
    phis = (+16.0, -16.0)  # two adjacent PMTs
    SPAN = 40

    # air interior backdrop
    gx, gy = glass_arc(INNER_GLASS, 0, span=SPAN)
    ax.fill(np.r_[gx, 120, 120], np.r_[gy, gy[-1], gy[0]], color=AIR, alpha=0.18, zorder=0)
    ax.plot([], [], color=AIR, lw=6, alpha=0.6, label='air (dom_inner)')

    # full band 187..203, then carve the domes out of it
    t = np.radians(np.linspace(-SPAN, SPAN, 500))
    out_arc = np.array([[INNER_GLASS * math.cos(a), INNER_GLASS * math.sin(a)] for a in t])
    in_arc  = np.array([[R_SHELL_IN * math.cos(a),  R_SHELL_IN * math.sin(a)]  for a in t])
    ax.fill(np.vstack([out_arc, in_arc[::-1]])[:, 0],
            np.vstack([out_arc, in_arc[::-1]])[:, 1],
            color=GOLD, alpha=0.6, zorder=1, label='gel = shell(187–203) − domes')

    for phi in phis:
        up, down = place_profile(outer_prof, MOUNT, phi)
        px = np.concatenate([up[:, 0], down[::-1, 0]])
        py = np.concatenate([up[:, 1], down[::-1, 1]])
        ax.fill(px, py, color='white', zorder=2)              # erase gel under the PMT
        ax.fill(px, py, color=GREEN, alpha=0.12, zorder=3)
        ax.plot(up[:, 0], up[:, 1], '-', color=GREEN, lw=2, zorder=4)
        ax.plot(down[:, 0], down[:, 1], '-', color=GREEN, lw=2, zorder=4)

    ax.plot(*glass_arc(INNER_GLASS, 0, span=SPAN), '-', color=BLUE, lw=1.8, zorder=5,
            label='inner glass (203)')
    ib = np.radians(np.linspace(-SPAN, SPAN, 300))
    ax.plot(R_SHELL_IN * np.cos(ib), R_SHELL_IN * np.sin(ib), ':', color='#7d6608',
            lw=1.6, zorder=5, label='shell inner bound (187)')

    aA = math.radians(16)
    ax.annotate('apex touches\nglass (203)', (203 * math.cos(aA), 203 * math.sin(aA)),
                xytext=(158, 70), fontsize=9, ha='center', arrowprops=dict(arrowstyle='->'))
    ax.annotate('thin meniscus gel\nwraps each dome', (192, 40),
                xytext=(142, 52), fontsize=9, color='#7d6608', ha='center',
                arrowprops=dict(arrowstyle='->', color=GOLD))
    ax.annotate('full gel in the valley\nbetween photocathodes', (197, 0),
                xytext=(150, 0), fontsize=9, color='#7d6608', ha='center',
                arrowprops=dict(arrowstyle='->', color=GOLD))
    ax.annotate('side walls = the dome faces\n= Boolean (shell − domes)',
                (186, -30), xytext=(131, -52), fontsize=9, color='#1b4f72', ha='left',
                arrowprops=dict(arrowstyle='->', color=BLUE))

    ax.set_title('Fat meniscus shell: gel fills 187→203 up to the dome faces', fontsize=12, pad=14)
    ax.set_xlabel('radial distance from DOM centre (mm)'); ax.set_ylabel('lateral (mm)')
    ax.set_aspect('equal'); ax.grid(alpha=0.25)
    ax.legend(loc='lower right', fontsize=8, framealpha=0.92)
    ax.set_xlim(128, 214); ax.set_ylim(-82, 82)

    plt.tight_layout(); plt.savefig(outfile, dpi=130, bbox_inches='tight'); plt.close()
    return outfile


if __name__ == "__main__":
    out = main()
    out2 = plot_fat_shell()
    print("wrote", out, "and", out2)
    print("knobs: GEL_THICK=%.1f  MOUNT_RECESSED=%.1f  half-angle=%.0f°"
          % (GEL_THICK, MOUNT_RECESSED, PHOTOCATHODE_HALF_ANGLE))
