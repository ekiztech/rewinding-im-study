# Compensating repair-induced core loss in rewound induction motors — code and data

Scripts and results for the paper *"Compensating repair-induced core loss in rewound induction motors by constrained
winding redesign: a loss-ratio criterion for repair shops"* (M. Ekiz, in preparation).

Reference machine: **CREATOR Case induction motor**, TU Graz (3.7 kW, 4 poles, 36/28 slots).
Machine data: K. Heidarikani, P. K. Dhakal, A. Mütze, R. Seebacher, *CREATOR Case: Induction Motor Data*,
TU Graz, 2024, doi:10.3217/kh3v7-dhn98 (CC BY-NC 4.0). Download it from the repository; it is not redistributed here
except for the small derived inputs listed below.

## Structure

| Folder | Content |
| --- | --- |
| `fea/` | FEMM 4.2 + pyFEMM model (run on Windows). Comments are in Turkish; file roles below. |
| `results/` | CSV outputs of the FEA runs used in the paper |
| `analysis/` | De-embedding, FEA-calibrated equivalent-circuit model, design scans, Monte Carlo, figures (any OS) |
| `figures/` | Generated figures (PNG + PDF) |

### `fea/` (run in this order; each writes a CSV next to itself)

| Script | Purpose | Output |
| --- | --- | --- |
| `im_geometry.py` | Parametric 2D geometry rebuilt from the DXF (one slot, rotated) | – |
| `adim2_dogrulama.py` | Model builder (materials, winding, circuits) + nominal validation | `sonuclar_adim2.csv` |
| `adim3_kalibrasyon.py`, `adim4_kalibrasyon2.py`, `adim4b_doluluk.py` | Calibration sweeps (air gap, B–H scaling, lamination fill, cage conductivity) | `sonuclar_adim3/4/4b.csv` |
| `adim5_son_dogrulama.py` | Calibrated model: no-load 50 Hz (11 points), locked rotor 10–50 Hz | `sonuclar_adim5.csv` |
| `adim6_adim_varyantlari.py` | Coil pitch 8/9 and 9/9: no-load curves | `sonuclar_adim6.csv` |
| `adim7_rotor_yeniden.py` | Locked rotor with corrected cage conductivity (29.3 MS/m), pitches 7/8/9 | `sonuclar_adim7.csv` |
| `adim8_harmonik_hasar.py` | Air-gap flux density samples (3 pitches) and degraded-permeability magnetizing curves | `sonuclar_adim8_bgap.csv`, `sonuclar_adim8_hasar.csv` |

Derived inputs: `bh_m800_50a.json` (M800-50A B–H from the dataset, extended above 1.7 T), `meas50.json`
(50 Hz no-load and locked-rotor points), `lr_freq.json` (locked-rotor points at 10–50 Hz).

Calibrated parameters: effective air gap 0.44 mm, effective lamination fill 0.87, cage conductivity 29.3 MS/m
(end-ring folded in, K_ring = 1.263), 3D end-effect inductance 2.91 mH (stator) + 3.15 mH (rotor).

### `analysis/` (Python 3, numpy, matplotlib; run from inside the folder)

| Script | Purpose |
| --- | --- |
| `deembed.py`, `deembed_run.py` | Separates the magnetizing branch from locked-rotor R and X (measured and FEA) → `deembedded2.json` |
| `hybrid2.py` | FEA-calibrated equivalent-circuit model (true parameters, three pitches, losses at 75 °C) |
| `scan2.py` | Discrete design scan: turns 14–20, 1–3 strands, IEC 60317 wires, 4 damage scenarios, 3 pitches |
| `compmap2.py` | Compensation map (Fig.) |
| `fig_pitch.py`, `harmloss.py` | Pitch comparison incl. belt-harmonic cage-loss estimate |
| `general.py`, `fig_general.py` | First version of the generality test (superseded by `rule2.py`/`rule3.py`, 12 loss splits) |
| `mc.py`, `fig_mc.py` | Monte Carlo uncertainty (N = 500, seed 2026; `mc.json` included) |
| `deemb_fig.py` | Error of skipping de-embedding and its effect on rated-point prediction |
| `perm_damage.py` | Core damage with reduced permeability + extended shop rule |
| `fig_noload.py`, `fig_lr2.py` | Validation figures |

Reproduce all analysis results and figures: `python deembed_run.py` first, then any other script.

## Licence

Code: MIT. Results and derived data: CC BY-NC 4.0 (follows the licence of the TU Graz dataset).

## Revision 1 (after internal review)

New analysis scripts (all in `analysis/`, run with Python 3, numpy, scipy, matplotlib):

| Script | What it does | Output |
| --- | --- | --- |
| `rule2.py` | Shop criterion, first-order and with flux-feedback term, against the full model on 12 synthetic loss splits | `rule2.json` |
| `ecgen.py` | IE2/IE3 motors (7.5–110 kW) built from published loss shares; criterion check, compensable damage, part load | `ecgen.json` |
| `refrev.py` (+ `refrev_lib.py`) | Reference motor: insulated lay factor, starting-current and MLT sensitivity, thermal feedback, part load, rated-current consistency | `refrev.json` |
| `mc2.py` | Monte Carlo with parameter, model-form and loop-test uncertainty (N = 500) | `mc2.json` |
| `figs_rev.py` | Figures: rule accuracy, pitch comparison, part load, Monte Carlo | `figures/` |

`fea9a.py` and `fea9b.py` analyse the FEA run 9 results (`results/sonuclar_adim9A.csv`, `results/sonuclar_adim9B2.csv`): leakage split, saturation at start, and harmonic cage losses.

New FEA script `fea/adim9_ek_fea.py`: saturated locked-rotor runs with an energy-based leakage split (part A) and rotating-rotor magnetostatic runs at rated load for harmonic cage losses (part B).

## Revision 2

| Script | What it does | Output |
| --- | --- | --- |
| `rule3.py` | Criterion with second-order terms (Eq. 6) against the full model, synthetic splits, IE motors, part load | `rule3.json` |
| `ecgen_sat.py` | IE2/IE3 motors with a saturating magnetizing branch | `ecgen_sat.json` |
| `dkmax_unc.py` | Absorbable core-loss increase with input uncertainty | `dkmax_unc.json` |
| `easa.py`, `fig_easa.py` | Measured before/after-rewind loss segregations of 20 motors (EASA/AEMT 2019 report, Tables 1-2 and 1-4; data in `data_ext/`) | `easa.json`, Fig. 9 |
| `scan_sat.py`, `fig_pitch_sat.py` | Design scan with saturated starting parameters; pitch comparison | `scan_sat.json`, Fig. 6 |
| `fig_rule.py` | Fig. 8 | |

Note: the FEA scripts in `fea/` carry Turkish comments; each file's header describes its purpose and outputs.
The calibration sweeps of steps 2–4 were exploratory; their final values are in Table 2 and reproduced by `adim5_son_dogrulama.py`.
