#!/usr/bin/env python3
"""
exotemp.py
==========

Surface temperatures and temperature maps of rocky exoplanets.

Physics
-------
* Insolation from stellar luminosity and orbital distance.
* Bond albedo A and heat-redistribution efficiency eps (0 = none, 1 = full),
  following Cowan & Agol (2011):
      T_day^4   = (S/sigma) (1 - A) (2/3 - 5 eps / 12)
      T_night^4 = (S/sigma) (1 - A) eps / 4
  eps = 1 recovers the familiar equilibrium temperature T_eq.
* A single-layer grey-atmosphere greenhouse model with infrared absorptivity g
  (0 = no atmosphere, 1 = opaque):  T_s = T_e (2 / (2 - g))^(1/4).
  Earth: T_e = 255 K and g = 0.78 give T_s = 288 K.
* Surface temperature maps of tidally locked planets: each point re-radiates
  the fraction (1 - eps) of the starlight it receives locally, and the fraction
  eps is spread evenly over the whole planet:
      sigma T_e^4(lat, lon) = S (1 - A) [ (1 - eps) max(cos Z, 0) + eps / 4 ]
  where Z is the angle from the substellar point. This conserves energy exactly.
* Derived products: the fraction of the surface able to hold liquid water,
  and the planet's thermal phase curve (day-night contrast).

Author: Javed Iqbal (https://github.com/Jiqbal-astro)
"""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

SIGMA = 5.670374419e-8          # Stefan-Boltzmann constant [W m^-2 K^-4]
SOLAR_CONSTANT = 1361.0         # Earth insolation [W m^-2] (Kopp & Lean 2011)
T_FREEZE, T_BOIL = 273.15, 373.15

DATA_DIR = Path(__file__).resolve().parent / "data"


# --------------------------------------------------------------------------
# Global (zero-dimensional) temperatures
# --------------------------------------------------------------------------
def insolation(luminosity, a_au):
    """Stellar flux at the planet [W m^-2] for luminosity in L_sun and a in AU."""
    return SOLAR_CONSTANT * np.asarray(luminosity, dtype=float) / np.asarray(a_au, dtype=float) ** 2


def equilibrium_temperature(flux, albedo):
    """Planet-averaged equilibrium temperature [K] (full heat redistribution)."""
    return (np.asarray(flux) * (1 - np.asarray(albedo)) / (4 * SIGMA)) ** 0.25


def dayside_temperature(flux, albedo, eps):
    """Dayside temperature [K] for redistribution efficiency eps (Cowan & Agol 2011)."""
    eps = np.clip(eps, 0, 1)
    return (np.asarray(flux) * (1 - np.asarray(albedo)) / SIGMA * (2 / 3 - 5 * eps / 12)) ** 0.25


def nightside_temperature(flux, albedo, eps):
    """Nightside temperature [K] for redistribution efficiency eps (Cowan & Agol 2011)."""
    eps = np.clip(eps, 0, 1)
    return (np.asarray(flux) * (1 - np.asarray(albedo)) / SIGMA * eps / 4) ** 0.25


def greenhouse_factor(g):
    """Surface warming factor T_s / T_e of a single-layer grey atmosphere."""
    g = np.asarray(g, dtype=float)
    if np.any((g < 0) | (g >= 2)):
        raise ValueError("Infrared absorptivity g must satisfy 0 <= g < 2 (physically 0-1).")
    return (2.0 / (2.0 - g)) ** 0.25


def surface_temperature(t_emission, g):
    """Surface temperature [K] below a single-layer atmosphere with absorptivity g."""
    return np.asarray(t_emission) * greenhouse_factor(g)


# --------------------------------------------------------------------------
# Surface maps
# --------------------------------------------------------------------------
def latlon_grid(n_lat=90, n_lon=180):
    """Cell-centred latitude/longitude grid [rad]; substellar point at (0, 0)."""
    lat = np.linspace(-np.pi / 2, np.pi / 2, n_lat + 1)
    lon = np.linspace(-np.pi, np.pi, n_lon + 1)
    lat_c = 0.5 * (lat[:-1] + lat[1:])
    lon_c = 0.5 * (lon[:-1] + lon[1:])
    LON, LAT = np.meshgrid(lon_c, lat_c)
    # exact area of each cell, normalised so the sphere sums to 1
    area = np.outer(np.diff(np.sin(lat)), np.diff(lon)) / (4 * np.pi)
    return LAT, LON, area


def temperature_map(flux, albedo, eps, g=0.0, grid=None):
    """Surface temperature map [K] of a tidally locked planet.

    Returns (T, LAT, LON, area) with the substellar point at lat = lon = 0.
    """
    LAT, LON, area = grid if grid is not None else latlon_grid()
    cos_z = np.cos(LAT) * np.cos(LON)
    absorbed = flux * (1 - albedo) * ((1 - eps) * np.clip(cos_z, 0, None) + eps / 4)
    t_e = (absorbed / SIGMA) ** 0.25
    return surface_temperature(t_e, g), LAT, LON, area


def liquid_water_fraction(T, area, t_min=T_FREEZE, t_max=T_BOIL):
    """Fraction of the planet's surface area with t_min <= T <= t_max."""
    return float(np.sum(area[(T >= t_min) & (T <= t_max)]))


def phase_curve(T, LAT, LON, area, phases):
    """Disk-integrated thermal emission versus orbital phase (edge-on orbit).

    Phase 0 = transit (nightside faces the observer), 0.5 = secondary eclipse
    (dayside faces the observer). Returned flux is normalised to its maximum.
    """
    phases = np.atleast_1d(phases)
    emission = SIGMA * T ** 4 * area
    flux = []
    for ph in phases:
        lon_obs = np.pi + 2 * np.pi * ph          # sub-observer longitude
        mu = np.cos(LAT) * np.cos(LON - lon_obs)
        flux.append(np.sum(emission * np.clip(mu, 0, None)))
    flux = np.asarray(flux)
    return flux / flux.max()


# --------------------------------------------------------------------------
# Data
# --------------------------------------------------------------------------
def load_planets(path=DATA_DIR / "planets.csv"):
    df = pd.read_csv(path)
    df["flux_w_m2"] = insolation(df.luminosity_lsun, df.semi_major_axis_au)
    df["flux_searth"] = df.flux_w_m2 / SOLAR_CONSTANT
    return df


def summary_table(planets, albedo=0.3, g=0.78):
    """Key temperatures for every planet, for the given albedo and greenhouse."""
    rows = []
    for _, p in planets.iterrows():
        S = p.flux_w_m2
        t_eq = float(equilibrium_temperature(S, albedo))
        row = {
            "planet": p.planet,
            "insolation_searth": p.flux_searth,
            "T_eq_K": t_eq,
            "T_surface_global_K": float(surface_temperature(t_eq, g)),
        }
        if p.tidally_locked:
            row["T_day_bare_rock_K"] = float(dayside_temperature(S, albedo, 0.0))
            T, LAT, LON, area = temperature_map(S, albedo, 0.3, g)
            row["T_substellar_eps0.3_K"] = float(T.max())
            row["water_fraction_eps0.3"] = liquid_water_fraction(T, area)
        rows.append(row)
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# Figures
# --------------------------------------------------------------------------
def plot_parameter_maps(planets, out, g=0.78):
    """Dayside surface temperature over Bond albedo and redistribution efficiency."""
    sel = planets[planets.tidally_locked].reset_index(drop=True)
    A = np.linspace(0.0, 0.6, 121)
    E = np.linspace(0.0, 1.0, 121)
    EE, AA = np.meshgrid(E, A)
    ncol = 4
    nrow = int(np.ceil(len(sel) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 3.6 * nrow), sharex=True, sharey=True)
    axes = np.atleast_1d(axes).ravel()
    levels = np.arange(150, 551, 25)
    for ax, (_, p) in zip(axes, sel.iterrows()):
        T = surface_temperature(dayside_temperature(p.flux_w_m2, AA, EE), g)
        cs = ax.contourf(EE, AA, np.clip(T, levels[0], levels[-1]), levels=levels, cmap="inferno")
        lines = ax.contour(EE, AA, T, levels=[T_FREEZE, T_BOIL], colors=["#9ecae1", "white"], linewidths=1.4)
        ax.clabel(lines, fmt={T_FREEZE: "273 K", T_BOIL: "373 K"}, fontsize=8)
        ax.set_title(f"{p.planet}  ({p.flux_searth:.2f} S⊕)", fontsize=10)
    for ax in axes[len(sel):]:
        ax.axis("off")
    for ax in axes[::ncol]:
        ax.set_ylabel("Bond albedo")
    for ax in axes[-ncol:]:
        ax.set_xlabel("Heat redistribution ε")
    fig.colorbar(cs, ax=axes.tolist(), label="Dayside surface temperature [K]", shrink=0.9)
    fig.suptitle(f"Dayside surface temperature (greenhouse absorptivity g = {g})", fontsize=13)
    fig.savefig(out, dpi=140, bbox_inches="tight")
    plt.close(fig)


def plot_surface_maps(planet, out, albedo=0.3, g=0.78, eps_values=(0.0, 0.3, 1.0)):
    """Mollweide maps of the surface temperature for several redistribution efficiencies."""
    fig = plt.figure(figsize=(5.2 * len(eps_values), 3.6))
    grid = latlon_grid(90, 180)
    maps = [temperature_map(planet.flux_w_m2, albedo, e, g, grid) for e in eps_values]
    vmax = max(m[0].max() for m in maps)
    for i, (eps, (T, LAT, LON, area)) in enumerate(zip(eps_values, maps), start=1):
        ax = fig.add_subplot(1, len(eps_values), i, projection="mollweide")
        mesh = ax.pcolormesh(LON, LAT, T, cmap="inferno", vmin=0, vmax=vmax, shading="auto")
        ax.contour(LON, LAT, T, levels=[T_FREEZE, T_BOIL], colors=["#9ecae1", "white"], linewidths=1.2)
        frac = liquid_water_fraction(T, area)
        ax.set_title(f"ε = {eps:g}   |   liquid-water area {frac:.0%}", fontsize=10)
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.grid(alpha=0.3)
    cb = fig.colorbar(mesh, ax=fig.axes, orientation="horizontal", fraction=0.06, pad=0.08, shrink=0.6)
    cb.set_label("Surface temperature [K]   (blue line 273 K, white line 373 K)")
    fig.suptitle(f"{planet.planet}: surface temperature maps (A = {albedo}, g = {g}); substellar point at centre",
                 fontsize=12)
    fig.savefig(out, dpi=140, bbox_inches="tight")
    plt.close(fig)


def plot_water_fraction(planets, out, g=0.78):
    """Liquid-water surface fraction versus redistribution efficiency."""
    sel = planets[planets.tidally_locked]
    eps = np.linspace(0, 1, 61)
    grid = latlon_grid(60, 120)
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.3), sharey=True)
    for ax, albedo in zip(axes, (0.1, 0.3, 0.5)):
        for i, (_, p) in enumerate(sel.iterrows()):
            frac = [liquid_water_fraction(temperature_map(p.flux_w_m2, albedo, e, g, grid)[0], grid[2]) for e in eps]
            ax.plot(eps, np.array(frac) * 100, lw=1.8, ls=("-", "--", ":")[i % 3], label=p.planet)
        ax.set_title(f"Bond albedo A = {albedo}")
        ax.set_xlabel("Heat redistribution ε")
        ax.grid(alpha=0.3)
    axes[0].set_ylabel("Surface area with 273–373 K [%]")
    axes[-1].legend(fontsize=8, loc="center left", bbox_to_anchor=(1.02, 0.5))
    fig.suptitle(f"How much of each planet could hold liquid water? (g = {g})", fontsize=13)
    fig.savefig(out, dpi=140, bbox_inches="tight")
    plt.close(fig)


def plot_phase_curves(planet, out, albedo=0.3, g=0.0, eps_values=(0.0, 0.1, 0.3, 0.6, 1.0)):
    """Thermal phase curves and day-night contrast."""
    phases = np.linspace(0, 1, 181)
    grid = latlon_grid(60, 120)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 4.3))
    cmap = plt.get_cmap("viridis")
    for i, e in enumerate(eps_values):
        T, LAT, LON, area = temperature_map(planet.flux_w_m2, albedo, e, g, grid)
        ax1.plot(phases, phase_curve(T, LAT, LON, area, phases), color=cmap(i / (len(eps_values) - 1)),
                 lw=2, label=f"ε = {e:g}")
    ax1.axvline(0.5, color="0.5", ls=":", lw=1)
    ax1.text(0.5, 0.02, " eclipse", color="0.4", fontsize=8)
    ax1.text(0.0, 0.02, " transit", color="0.4", fontsize=8)
    ax1.set_xlabel("Orbital phase")
    ax1.set_ylabel("Thermal emission (normalised)")
    ax1.set_title(f"{planet.planet}: thermal phase curves")
    ax1.legend(fontsize=8, loc="upper left")
    ax1.grid(alpha=0.3)

    eps = np.linspace(0, 1, 101)
    td = dayside_temperature(planet.flux_w_m2, albedo, eps)
    tn = nightside_temperature(planet.flux_w_m2, albedo, eps)
    ax2.plot(eps, td, lw=2, color="tab:red", label="Dayside")
    ax2.plot(eps, tn, lw=2, color="tab:blue", label="Nightside")
    ax2.axhline(float(equilibrium_temperature(planet.flux_w_m2, albedo)), color="k", ls="--", lw=1, label="T_eq")
    ax2.set_xlabel("Heat redistribution ε")
    ax2.set_ylabel("Brightness temperature [K]")
    ax2.set_title("Day–night contrast (Cowan & Agol 2011)")
    ax2.legend(fontsize=8)
    ax2.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(out, dpi=140, bbox_inches="tight")
    plt.close(fig)


def plot_3d_surface(planet, out, g=0.78):
    """3D view of dayside surface temperature over albedo and redistribution."""
    A = np.linspace(0.0, 0.6, 80)
    E = np.linspace(0.0, 1.0, 80)
    EE, AA = np.meshgrid(E, A)
    T = surface_temperature(dayside_temperature(planet.flux_w_m2, AA, EE), g)
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection="3d")
    surf = ax.plot_surface(EE, AA, T, cmap="inferno", edgecolor="none")
    ax.set_xlabel("Heat redistribution ε")
    ax.set_ylabel("Bond albedo")
    ax.set_zlabel("T [K]")
    ax.set_title(f"{planet.planet}: dayside surface temperature (g = {g})")
    fig.colorbar(surf, shrink=0.6, label="Temperature [K]")
    fig.savefig(out, dpi=140, bbox_inches="tight")
    plt.close(fig)


# --------------------------------------------------------------------------
# Command line
# --------------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser(description="Surface temperatures and temperature maps of rocky exoplanets.")
    p.add_argument("-o", "--output", default="figures", help="output directory (default: figures)")
    p.add_argument("--planets", default=DATA_DIR / "planets.csv", help="planet catalogue CSV")
    p.add_argument("--albedo", type=float, default=0.3, help="Bond albedo for maps and table (default 0.3)")
    p.add_argument("--greenhouse", type=float, default=0.78,
                   help="infrared absorptivity g of the atmosphere, 0-1 (default 0.78, Earth-like)")
    p.add_argument("--map-planet", default="Proxima Cen b", help="planet for maps and phase curves")
    args = p.parse_args()

    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    planets = load_planets(args.planets)
    table = summary_table(planets, args.albedo, args.greenhouse)
    table.to_csv(out / "temperature_summary.csv", index=False, float_format="%.4g")
    pd.set_option("display.width", 160)
    print(table.to_string(index=False, float_format=lambda x: f"{x:.3g}"))

    target = planets[planets.planet == args.map_planet]
    if target.empty:
        raise SystemExit(f"Planet '{args.map_planet}' not in catalogue")
    target = target.iloc[0]
    plot_parameter_maps(planets, out / "dayside_parameter_maps.png", args.greenhouse)
    plot_surface_maps(target, out / "surface_maps.png", args.albedo, args.greenhouse)
    plot_water_fraction(planets, out / "liquid_water_fraction.png", args.greenhouse)
    plot_phase_curves(target, out / "phase_curves.png", args.albedo)
    plot_3d_surface(target, out / "dayside_3d.png", args.greenhouse)
    print(f"\nFigures and table saved to {out}/")


if __name__ == "__main__":
    main()
