"""Physics checks against textbook and published values. Run with: pytest"""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import exotemp as et  # noqa: E402


def test_earth_equilibrium_and_surface_temperature():
    # 255 K without greenhouse; single-layer model with g = 0.78 gives ~288 K
    t_eq = et.equilibrium_temperature(et.SOLAR_CONSTANT, 0.3)
    assert t_eq == pytest.approx(255, abs=1)
    assert et.surface_temperature(t_eq, 0.78) == pytest.approx(288, abs=1.5)


def test_greenhouse_always_warms():
    assert et.greenhouse_factor(0.0) == 1.0
    g = np.linspace(0, 1, 11)
    assert np.all(np.diff(et.greenhouse_factor(g)) > 0)
    with pytest.raises(ValueError):
        et.greenhouse_factor(-0.1)


def test_full_redistribution_recovers_equilibrium_temperature():
    S, A = 900.0, 0.3
    t_eq = et.equilibrium_temperature(S, A)
    assert et.dayside_temperature(S, A, 1.0) == pytest.approx(t_eq)
    assert et.nightside_temperature(S, A, 1.0) == pytest.approx(t_eq)


def test_less_redistribution_means_hotter_dayside():
    eps = np.linspace(0, 1, 21)
    td = et.dayside_temperature(900, 0.3, eps)
    tn = et.nightside_temperature(900, 0.3, eps)
    assert np.all(np.diff(td) < 0) and np.all(np.diff(tn) > 0)
    # bare rock (eps = 0): T_day = (8/3)^(1/4) T_eq
    assert td[0] / et.equilibrium_temperature(900, 0.3) == pytest.approx((8 / 3) ** 0.25)


@pytest.mark.parametrize("eps", [0.0, 0.25, 0.5, 1.0])
def test_map_conserves_energy(eps):
    S, A = 1200.0, 0.2
    T, LAT, LON, area = et.temperature_map(S, A, eps, g=0.0, grid=et.latlon_grid(180, 360))
    emitted = np.sum(et.SIGMA * T**4 * area) * 4      # per unit cross-section, x4 for sphere/disk
    assert emitted == pytest.approx(S * (1 - A), rel=2e-3)


@pytest.mark.parametrize("eps", [0.0, 0.3, 0.7, 1.0])
def test_map_eclipse_brightness_matches_cowan_agol(eps):
    # Disk-integrated dayside emission of the map equals the analytic
    # Cowan & Agol (2011) dayside temperature for every eps.
    S, A = 1000.0, 0.3
    T, LAT, LON, area = et.temperature_map(S, A, eps, grid=et.latlon_grid(180, 360))
    mu = np.clip(np.cos(LAT) * np.cos(LON), 0, None)
    flux = np.sum(et.SIGMA * T**4 * mu * area) * 4 * np.pi   # integral over solid angle
    t_day = (flux / (np.pi * et.SIGMA)) ** 0.25
    assert t_day == pytest.approx(et.dayside_temperature(S, A, eps), rel=2e-3)


def test_grid_area_sums_to_one():
    assert et.latlon_grid(30, 60)[2].sum() == pytest.approx(1.0)


def test_phase_curve_flat_for_full_redistribution_and_peaks_at_eclipse():
    phases = np.linspace(0, 1, 41)
    grid = et.latlon_grid(60, 120)
    flat = et.phase_curve(et.temperature_map(900, 0.3, 1.0, grid=grid)[0], *grid, phases)
    assert np.ptp(flat) < 1e-3
    locked = et.phase_curve(et.temperature_map(900, 0.3, 0.0, grid=grid)[0], *grid, phases)
    assert phases[np.argmax(locked)] == pytest.approx(0.5)
    assert locked[0] < 0.01        # nightside at transit is dark


def test_liquid_water_fraction_limits():
    grid = et.latlon_grid(60, 120)
    hot = et.temperature_map(20 * et.SOLAR_CONSTANT, 0.3, 1.0, grid=grid)[0]
    assert et.liquid_water_fraction(hot, grid[2]) == 0.0
    earthlike = et.temperature_map(et.SOLAR_CONSTANT, 0.3, 1.0, g=0.78, grid=grid)[0]
    assert et.liquid_water_fraction(earthlike, grid[2]) == pytest.approx(1.0)


def test_catalogue_insolations_match_literature():
    df = et.load_planets().set_index("planet")
    assert df.loc["Earth", "flux_searth"] == pytest.approx(1.0)
    assert df.loc["Proxima Cen b", "flux_searth"] == pytest.approx(0.65, abs=0.02)   # Anglada-Escudé+ 2016
    assert df.loc["TRAPPIST-1 e", "flux_searth"] == pytest.approx(0.646, abs=0.01)   # Agol+ 2021
    assert df.loc["TOI-700 d", "flux_searth"] == pytest.approx(0.87, abs=0.03)       # Gilbert+ 2020


def test_trappist1e_equilibrium_temperature_matches_agol2021():
    # Agol et al. (2021) give T_eq = 249.7 K for zero albedo
    df = et.load_planets().set_index("planet")
    assert et.equilibrium_temperature(df.loc["TRAPPIST-1 e", "flux_w_m2"], 0.0) == pytest.approx(249.7, abs=2)
