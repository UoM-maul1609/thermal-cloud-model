#!/usr/bin/env python3
"""Run BAM activation sweeps for VOCALS aerosol and added NaCl spray.

Two experiments are performed:
  1. VOCALS aerosol only: w = 1e-2 ... 10 m s-1 (100 log-spaced points).
  2. VOCALS + NaCl spray: w = 0.3 m s-1 and q_NaCl = 1e-14 ... 1e-4
     kg kg-1 (50 log-spaced points).

Schemes compared:
  ARG, Ghosh-modified ARG, FN, FN + giant-CCN correction, and FN quadrature.

For the Ghosh scheme, GHOSH_SIGMA_MODE selects either the published-style
fixed accumulation-mode width or an experimental effective width diagnosed by
BAM from the combined aerosol PSD.

Set SV_FLAG = 1 below to include the 10-bin semi-volatile organic distribution.
Set SV_FLAG = 0 to reproduce the non-semi-volatile sweeps.

BAM's public namelist interface expects aerosol number mixing ratios in
kg-1 dry air, semi-volatile organic loadings in microgram kg-1 dry air,
and returns activated droplet number in kg-1 dry air.
"""

from dataclasses import dataclass
from pathlib import Path
import getpass
import subprocess
import numpy as np
import matplotlib.pyplot as plt


# -----------------------------------------------------------------------------
# User settings
# -----------------------------------------------------------------------------
# This script is intended to live in:
#     <project_root>/python/
# and to be launched from the project root with:
#     python3 python/bam_vocals_nacl_sweeps_native_perkg.py
#
# Paths are based on the script location rather than the shell's current
# working directory, so the BAM executable is found reliably.
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent

BAM_EXE = PROJECT_ROOT / "main.exe"

# Put all temporary namelists, CSV files and figures outside the source tree:
#     /tmp/<username>/bam_sweep_output/
USER = getpass.getuser()
OUTPUT_DIR = Path("/tmp") / USER / "bam_sweep_output"

SHOW_PLOTS = True

# Cloud-base conditions used by the supplied BAM namelist.
P_TEST = 99344.0       # Pa
T_TEST = 280.15        # K

# ARG fitting coefficients used by the supplied BAM namelist.
A_EQ_7 = 0.21
B_EQ_7 = 1.58

# Ghosh et al. (2025) modified-ARG settings.
#
# GHOSH_SIGMA_MODE = 0:
#   Published-style option. Use one fixed accumulation-mode geometric standard
#   deviation sigma_g for all aerosol modes. For the VOCALS 154 nm mode,
#   sigma_g = exp(0.465) ~= 1.592.
#
# GHOSH_SIGMA_MODE = 1:
#   Experimental MCB extension. BAM diagnoses a single effective sigma_g from
#   the *combined* current PSD within GHOSH_DMIN ... GHOSH_DMAX and applies the
#   resulting Ghosh f, g and p coefficients to every mode.
GHOSH_SIGMA_MODE = 0
GHOSH_SIGMA_ACC = float(np.exp(0.465))  # sigma_g, not ln(sigma_g)
GHOSH_DMIN = 80.0e-9                    # m
GHOSH_DMAX = 1.0e-6                     # m

# Semi-volatile organic treatment.
# 0 = off, 1 = on.
SV_FLAG = 0

# Optional multiplier applied to every volatility-bin loading.  Leave at 1.0
# to use the values from the supplied BAM namelist.
SV_SCALE = 1.0

# 10-bin VBS used by the supplied BAM namelist.
# org_content is microgram organic / kg dry air.
SV_ORG_CONTENT = SV_SCALE * np.array(
    [0.005, 0.01, 0.02, 0.03, 0.06, 0.08, 0.16, 0.30, 0.42, 0.80],
    dtype=float,
)
SV_MOLW = np.full(10, 200.0e-3)          # kg mol-1
SV_RHO = np.full(10, 1500.0)             # kg m-3
SV_DELTA_H_VAP = np.full(10, 150.0)      # kJ mol-1
SV_NU = np.ones(10)
SV_LOG_C_STAR = np.arange(-6.0, 4.0, 1.0)  # log10(C* / (microgram m-3))

# Sweep definitions.
UPDRAFTS = np.logspace(-2.0, 1.0, 100)   # m s-1
SALT_SWEEP_W = 0.3                        # m s-1
NACL_MR = np.logspace(-14.0, -4.0, 50)   # kg NaCl / kg dry air

# Salt-spray PSD selector.
COOPER = 1
HARRISON_EFFERVESCENCE = 2
HARRISON_DELAVAL = 3
EDMUND_SH = 4
RAYLEIGH = 5
TAYLOR_CONE = 6
TOPAZ = 7
POWERCLOUD = 8

spray_method = EDMUND_SH


# -----------------------------------------------------------------------------
# Aerosol definitions
# -----------------------------------------------------------------------------
# VOCALS accumulation/Aitken PSD.  Number is supplied directly in BAM's
# native number-mixing-ratio unit, kg-1 dry air. sig is ln(sigma_g).
VOCALS_N = np.array([46.6469e6, 153.42e6, 166.77e6], dtype=float)
VOCALS_D = np.array([18e-9, 39e-9, 154e-9], dtype=float)            # m
VOCALS_LOGSIG = np.array([0.348, 0.354, 0.465], dtype=float)

# Existing BAM core chemistry for the VOCALS modes: ammonium-sulfate-like.
VOCALS_MOLW = np.full(3, 132.14e-3)   # kg mol-1
VOCALS_RHO = np.full(3, 1770.0)       # kg m-3
VOCALS_NU = np.full(3, 3.0)

# NaCl properties for every appended spray mode.
NACL_MOLW = 58.44e-3                  # kg mol-1
NACL_RHO = 2165.0                     # kg m-3
NACL_NU = 2.0


@dataclass(frozen=True)
class Scheme:
    label: str
    method_flag: int
    giant_flag: int


if GHOSH_SIGMA_MODE == 0:
    GHOSH_LABEL = "Ghosh ARG (fixed)"
elif GHOSH_SIGMA_MODE == 1:
    GHOSH_LABEL = "Ghosh ARG (effective)"
else:
    raise ValueError("GHOSH_SIGMA_MODE must be 0 (fixed) or 1 (effective)")

SCHEMES = (
    Scheme("ARG", 1, 0),
    Scheme(GHOSH_LABEL, 4, 0),
    Scheme("FN", 2, 0),
    Scheme("FN + giant", 2, 1),
    Scheme("FN quadrature", 3, 0),
)

_GHOSH_RANGE_WARNING_SHOWN = False


def get_spray_psd(method):
    """Return (relative number weights, ln(sigma_g), Dg [m], name)."""
    if method == COOPER:
        n = np.array([1.55, 194.0, 17.7]) * 1e6
        logsig = np.array([0.129, 0.625, 0.666])
        dm = np.array([0.0263, 0.0588, 0.269]) * 1e-6
        name = "COOPER"

    elif method == HARRISON_EFFERVESCENCE:
        n = np.array([94800.0, 245000.0, 137000.0]) * 1e6
        logsig = np.array([0.252, 0.497, 0.883])
        dm = np.array([0.0264, 0.0434, 0.0574]) * 1e-6
        name = "HARRISON_EFFERVESCENCE"

    elif method == HARRISON_DELAVAL:
        n = np.array([214000.0, 66800.0, 4800.0]) * 1e6
        logsig = np.array([0.703, 0.534, 0.900])
        dm = np.array([0.0156, 0.123, 0.398]) * 1e-6
        name = "HARRISON_DELAVAL"

    elif method == EDMUND_SH:
        # New version including the small-particle mode.
        n = np.array([250000.0, 789000.0, 2050000.0]) * 1e6
        logsig = np.array([0.200, 0.300, 0.751])
        dm = np.array([0.00941, 0.0168, 0.0467]) * 1e-6
        name = "EDMUND_SH"

    elif method == RAYLEIGH:
        n = np.array([1.0]) * 1e6
        logsig = np.array([0.25])
        dm = np.array([0.162]) * 1e-6
        name = "RAYLEIGH"

    elif method == TAYLOR_CONE:
        n = np.array([1.0]) * 1e6
        logsig = np.log(np.array([1.19]))
        dm = np.array([0.0826]) * 1e-6
        name = "TAYLOR_CONE"

    elif method == TOPAZ:
        n = np.array([3530.0, 190.0, 0.00378]) * 1e6
        logsig = np.array([0.424, 0.374, 0.010])
        dm = np.array([0.0848, 0.198, 0.915]) * 1e-6
        name = "TOPAZ"

    elif method == POWERCLOUD:
        n = np.array([11.59948504, 124.6248093, 51.36088692]) * 1e6
        logsig = np.log(np.array([1.163293663, 1.63031549, 1.693508788]))
        dm = np.array([15.76237815, 24.49137812, 139.9466964]) * 1e-9
        name = "POWERCLOUD"

    else:
        raise ValueError(f"Unknown spray_method={method}")

    # Only the relative mode weights are used when scaling to a specified
    # NaCl mass mixing ratio.
    weights = n / np.sum(n)
    return weights, logsig, dm, name


def nacl_number_from_mixing_ratio(q_nacl, weights, logsig, dm):
    """Scale a multimodal NaCl PSD to q_nacl [kg kg-1 dry air].

    For a lognormal mode with median diameter Dg and logsig=ln(sigma_g),
    <D^3> = Dg^3 exp(4.5 logsig^2).

    BAM's namelist number unit is kg-1 dry air. Therefore a NaCl mass
    mixing ratio q_nacl [kg kg-1] divided by the shape-mean particle mass
    directly gives total particle number [kg-1].
    """
    mean_particle_mass = (
        np.pi / 6.0
        * NACL_RHO
        * dm**3
        * np.exp(4.5 * logsig**2)
    )
    mean_mass_for_shape = np.sum(weights * mean_particle_mass)
    n_total = q_nacl / mean_mass_for_shape
    return n_total * weights


def append_nacl_to_vocals(q_nacl, weights, salt_logsig, salt_dm):
    """Return BAM modal arrays with the NaCl modes appended after VOCALS."""
    salt_n = nacl_number_from_mixing_ratio(q_nacl, weights, salt_logsig, salt_dm)
    nsalt = len(salt_n)

    n = np.concatenate((VOCALS_N, salt_n))
    d = np.concatenate((VOCALS_D, salt_dm))
    logsig = np.concatenate((VOCALS_LOGSIG, salt_logsig))
    molw = np.concatenate((VOCALS_MOLW, np.full(nsalt, NACL_MOLW)))
    rho = np.concatenate((VOCALS_RHO, np.full(nsalt, NACL_RHO)))
    nu = np.concatenate((VOCALS_NU, np.full(nsalt, NACL_NU)))
    return n, d, logsig, molw, rho, nu


def append_array(lines, name, values):
    """Write one namelist element per record for readable generated cases."""
    for i, value in enumerate(np.asarray(values), start=1):
        lines.append(f"  {name}({i}) = {float(value):.16e},")


def write_bam_namelist(path, scheme, w, n, d, logsig, molw, rho, nu):
    """Write a BAM namelist for one deterministic activation calculation."""
    nmode = len(n)

    if SV_FLAG not in (0, 1):
        raise ValueError("SV_FLAG must be 0 or 1")

    if SV_FLAG == 1:
        nsv = len(SV_ORG_CONTENT)
        sv_arrays = (
            SV_ORG_CONTENT,
            SV_MOLW,
            SV_RHO,
            SV_DELTA_H_VAP,
            SV_NU,
            SV_LOG_C_STAR,
        )
        if any(len(a) != nsv for a in sv_arrays):
            raise ValueError("All semi-volatile arrays must have the same length")
        if np.any(SV_ORG_CONTENT < 0.0):
            raise ValueError("SV_ORG_CONTENT must be non-negative")
    else:
        # BAM currently requires n_sv > 0 even when sv_flag=0, so retain one
        # harmless dummy bin in the disabled case.
        nsv = 1

    lines = [
        "&bulk_aerosol_setup",
        f"  n_mode = {nmode},",
        f"  n_sv = {nsv},",
        f"  sv_flag = {SV_FLAG},",
        f"  method_flag = {scheme.method_flag},",
        f"  giant_flag = {scheme.giant_flag},",
        f"  a_eq_7 = {A_EQ_7:.16e},",
        f"  b_eq_7 = {B_EQ_7:.16e},",
        f"  ghosh_sigma_mode = {GHOSH_SIGMA_MODE},",
        f"  ghosh_sigma_acc = {GHOSH_SIGMA_ACC:.16e},",
        f"  ghosh_dmin = {GHOSH_DMIN:.16e},",
        f"  ghosh_dmax = {GHOSH_DMAX:.16e},",
        "/",
        "&bulk_aerosol_spec",
    ]

    append_array(lines, "n_aer1", n)
    append_array(lines, "d_aer1", d)
    append_array(lines, "sig_aer1", logsig)
    append_array(lines, "molw_core1", molw)
    append_array(lines, "density_core1", rho)
    append_array(lines, "nu_core1", nu)

    if SV_FLAG == 1:
        append_array(lines, "org_content1", SV_ORG_CONTENT)
        append_array(lines, "molw_org1", SV_MOLW)
        append_array(lines, "density_org1", SV_RHO)
        append_array(lines, "delta_h_vap1", SV_DELTA_H_VAP)
        append_array(lines, "nu_org1", SV_NU)
        append_array(lines, "log_c_star1", SV_LOG_C_STAR)
    else:
        # Dummy values are ignored when sv_flag=0.
        lines.extend(
            [
                "  org_content1(1) = 0.0,",
                "  molw_org1(1) = 2.0e-1,",
                "  density_org1(1) = 1500.0,",
                "  delta_h_vap1(1) = 150.0,",
                "  nu_org1(1) = 1.0,",
                "  log_c_star1(1) = 0.0,",
            ]
        )

    lines.extend(
        [
            f"  p_test = {P_TEST:.16e},",
            f"  t_test = {T_TEST:.16e},",
            f"  w_test = {float(w):.16e},",
            "  rand_dist = .false.,",
            "  n_rand = 1,",
            "  mean_w = 0.0,",
            "  sigma_w = 1.0,",
            "/",
        ]
    )

    path.write_text("\n".join(lines) + "\n")

def parse_ndrop(stdout, nmode):
    """Extract sum(act_frac1*n_aer1) from BAM's final print record."""
    for line in reversed(stdout.splitlines()):
        values = []
        for token in line.replace("D", "E").replace("d", "e").split():
            try:
                values.append(float(token))
            except ValueError:
                pass

        # final BAM record is:
        # w, act_frac(1:nmode), total_activated_number, total_fraction
        if len(values) >= nmode + 3:
            return values[-2]

    raise RuntimeError("Could not find activated drop number in BAM output:\n" + stdout)


def run_bam(case_file, scheme, w, n, d, logsig, molw, rho, nu):
    global _GHOSH_RANGE_WARNING_SHOWN

    write_bam_namelist(case_file, scheme, w, n, d, logsig, molw, rho, nu)

    proc = subprocess.run(
        [str(BAM_EXE), str(case_file)],
        cwd=str(BAM_EXE.parent),
        text=True,
        capture_output=True,
    )

    if (
        "BAM warning: diagnosed Ghosh sigma_g=" in proc.stdout
        and not _GHOSH_RANGE_WARNING_SHOWN
    ):
        print(
            "  WARNING: the experimental Ghosh effective width falls outside "
            "the published 1.4-2.1 calibration range for at least one case."
        )
        _GHOSH_RANGE_WARNING_SHOWN = True

    if proc.returncode != 0:
        raise RuntimeError(
            f"BAM failed for {scheme.label}, w={w:g} m/s\n"
            f"stdout:\n{proc.stdout}\n"
            f"stderr:\n{proc.stderr}\n"
            "BAM returned a non-zero exit status. Check the generated namelist "
            "and BAM diagnostics above."
        )

    return parse_ndrop(proc.stdout, len(n))


def ghosh_suffix():
    """Keep the default fixed-Ghosh filenames unchanged; tag effective runs."""
    return "_ghosh_effective" if GHOSH_SIGMA_MODE == 1 else ""


def output_suffix():
    """Filename suffix for the selected Ghosh and semi-volatile options."""
    return ghosh_suffix() + ("_sv" if SV_FLAG == 1 else "")


def sv_title():
    return " + semi-volatiles" if SV_FLAG == 1 else ""


def save_table(path, xname, x, results):
    names = [xname] + [scheme.label.replace(" ", "_") + "_kg-1" for scheme in SCHEMES]
    columns = [x] + [results[scheme.label] for scheme in SCHEMES]
    np.savetxt(path, np.column_stack(columns), delimiter=",", header=",".join(names), comments="")


def updraft_sweep(case_file):
    print("Running VOCALS updraft sweep...")
    results = {}

    for scheme in SCHEMES:
        print(f"  {scheme.label}")
        ndrop = np.empty_like(UPDRAFTS)
        for j, w in enumerate(UPDRAFTS):
            ndrop[j] = run_bam(
                case_file,
                scheme,
                w,
                VOCALS_N,
                VOCALS_D,
                VOCALS_LOGSIG,
                VOCALS_MOLW,
                VOCALS_RHO,
                VOCALS_NU,
            )
        results[scheme.label] = ndrop

    save_table(OUTPUT_DIR / f"vocals_updraft_sweep{output_suffix()}.csv", "w_m_s-1", UPDRAFTS, results)

    fig, ax = plt.subplots(figsize=(7.0, 5.0))
    for scheme in SCHEMES:
        ax.plot(UPDRAFTS, results[scheme.label], label=scheme.label)
    ax.set_xscale("log")
    ax.set_xlabel(r"Updraft speed (m s$^{-1}$)")
    ax.set_ylabel(r"Activated drop number (kg$^{-1}$ dry air)")
    ax.set_title(f"VOCALS aerosol activation{sv_title()}")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / f"vocals_updraft_sweep{output_suffix()}.png", dpi=200)
    return fig, results


def salt_sweep(case_file):
    weights, salt_logsig, salt_dm, spray_name = get_spray_psd(spray_method)
    print(f"Running VOCALS + NaCl sweep using {spray_name} PSD at w={SALT_SWEEP_W:g} m/s...")
    results = {}

    for scheme in SCHEMES:
        print(f"  {scheme.label}")
        ndrop = np.empty_like(NACL_MR)
        for j, q_nacl in enumerate(NACL_MR):
            n, d, logsig, molw, rho, nu = append_nacl_to_vocals(
                q_nacl, weights, salt_logsig, salt_dm
            )
            ndrop[j] = run_bam(
                case_file, scheme, SALT_SWEEP_W, n, d, logsig, molw, rho, nu
            )
        results[scheme.label] = ndrop

    stem = spray_name.lower()
    save_table(
        OUTPUT_DIR / f"vocals_nacl_{stem}_sweep{output_suffix()}.csv",
        "NaCl_kg_kg-1",
        NACL_MR,
        results,
    )

    fig, ax = plt.subplots(figsize=(7.0, 5.0))
    for scheme in SCHEMES:
        ax.plot(NACL_MR, results[scheme.label], label=scheme.label)
    ax.set_xscale("log")
    ax.set_xlabel(r"NaCl mixing ratio (kg kg$^{-1}$)")
    ax.set_ylabel(r"Activated drop number (kg$^{-1}$ dry air)")
    ax.set_title(f"VOCALS + NaCl ({spray_name}){sv_title()}, w = {SALT_SWEEP_W:g} m s$^{{-1}}$")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / f"vocals_nacl_{stem}_sweep{output_suffix()}.png", dpi=200)
    return fig, results


def main():
    if not BAM_EXE.exists():
        raise FileNotFoundError(
            f"Cannot find BAM executable: {BAM_EXE}\n"
            "Set BAM_EXE near the top of this script."
        )

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    if GHOSH_SIGMA_MODE == 0:
        print(
            f"Ghosh ARG: fixed sigma_acc = {GHOSH_SIGMA_ACC:.6g} "
            "(geometric standard deviation)"
        )
    else:
        print(
            f"Ghosh ARG: experimental effective sigma over "
            f"{GHOSH_DMIN*1e9:.3g}-{GHOSH_DMAX*1e9:.3g} nm"
        )

    if SV_FLAG == 1:
        print(
            f"Semi-volatiles ON: {len(SV_ORG_CONTENT)} bins, "
            f"total organic loading = {np.sum(SV_ORG_CONTENT):.6g} microgram kg-1"
        )
    else:
        print("Semi-volatiles OFF")

    case_file = OUTPUT_DIR / "_bam_sweep_case.in"

    fig1, _ = updraft_sweep(case_file)
    fig2, _ = salt_sweep(case_file)

    print(f"Results written to: {OUTPUT_DIR}")
    if SHOW_PLOTS:
        plt.show()
    else:
        plt.close(fig1)
        plt.close(fig2)


if __name__ == "__main__":
    main()
