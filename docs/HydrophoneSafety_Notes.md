# Hydrophone Safety Measurements — Context & Findings

_Last updated: 2026-08-13_

Working notes for hydrophone-based acoustic-output safety measurements (MI, I_sppa,
I_spta) of custom Verasonics ultrasound sequences, captured with an NI oscilloscope.

---

## 1. Setup

- **Goal:** measure MI and intensities for custom sequences (focused / wide-beam /
  diverging; ARF pushes for shear-wave elastography) and confirm they agree with the
  Verasonics-reported values.
- **DAQ:** NI **PCI-5112** digitizer — 8-bit, 100 MHz bandwidth, **100 MS/s real-time
  max**, 16–32 MB/ch. Driven by the **Scope Soft Front Panel (Scope-SFP 3.8)**, which
  writes `.hws` files.
- **Hydrophone chain:** Onda **HGL-0400** (SN 1037) + **AH-2010-025** preamp (SN 1109),
  calibrated 2024-06-10 (cal files in the repo `calibration/` folder; also baked into
  `HydrophoneTables.mat`).
- **Probes:** L11-5v (used for the MI validation below), **S5-1** (next).
- **MI pipeline of record:** `CalculateSafety.m` in
  `D:\Luuk van Knippenberg\Github\SWI\Mechanical index`.

---

## 2. `.hws` file format (NI Hierarchical Waveform Storage)

`.hws` **is HDF5** (magic bytes `‰HDF`). Read directly in MATLAB with `h5read` /
`h5readatt` — no NI drivers needed. Layout:

| Item | HDF5 path |
|---|---|
| Raw samples (int8 ADC codes) | `/wfm_group0/vectors/vector0/data` |
| Y scaling `[c0 c1]` → `V = c0 + c1·raw` | `/wfm_group0/axes/axis1/scale_coef` |
| Time axis `start`, `increment` (implicit) | attrs on `/wfm_group0/axes/axis0` |
| Full scope config (chan/hor/trig) | `/cfg_scope0/data/...` |
| Digitizer model | `/cfg_scope0/data/glb/public/names` |
| User notes | attr `note`/`notes` on `/wfm_group0/id` (+ mirrors) |

**Samples are raw integer codes, not volts** — the `scale_coef` step is essential.

---

## 3. Tools built (`D:\Luuk van Knippenberg\Claude\MI estimation\`)

| File | Purpose |
|---|---|
| `readHWS.m` | Read waveform `(t, y)` + `cfg` (sample rate, record length, n_records, bandwidth, impedance, coupling, digitizer, **notes**, **is_RIS**). Warns on RIS and on record-shorter-than-pulse. |
| `scanPlan.m` | Non-uniform 3-orthogonal-line scan sized to the beam (λ, lateral FWHM, DOF). Tiers via `anchorOnly` / `lines` subset. `center` = observed peak. Returns per-point trigger delay `tof_s`. |
| `measurementPlan.m` | Per-beam-type protocol: coarse **search** → recenter on peak → `.fine(peakXYZ)` (function handle) → 1-D voltage/pulse-length **sweeps** at the peak point. |
| `MI L11-5\PlotMI_L11_5.m` | MI-vs-TX-voltage analysis: calibration, impedance load-correction, derating, saturation-aware extrapolation. Options: `SpatialAvgCorr`/`Fnum`, `fMI`, `Method` (`singlefreq` default / `broadband`). |
| `MI L11-5\spatialAvgFactor.m` | Finite-aperture (spatial-averaging) correction, jinc focal-plane model. |
| `MI L11-5\voltageToPressureBroadband.m` | Frequency-domain deconvolution to pressure (hydrophone magnitude + **amplifier phase** + capacitive divider). |

**Reproducibility:** `PlotMI_L11_5` defaults (`Method="singlefreq"`, `CF=1`,
`fMI=measured`) reproduce the original single-frequency analysis exactly; `broadband`
writes a separate `_broadband.png` and never overwrites the canonical figure.

---

## 4. Measurement strategy

- **Measure at the OBSERVED peak, not the geometric focus.** The pressure max is
  displaced from the nominal focus by diffraction (focal shift, proximal), nonlinearity
  (further proximal), and aberration/lens (lateral + elevation). Diverging/wide beams
  have no focal peak — theirs sits in the near field. Workflow: coarse **search** →
  recenter fine scan on the found peak. Track **peak-negative pressure** for MI,
  **pulse-intensity-integral** for I_sppa.
- **Split parameters:**
  - _Changes field shape_ (beam type, #elements/aperture, focal depth) → needs a
    measurement (validate the simulation once per geometry).
  - _Scales predictably_ (transmit voltage, pulse length, PRF) → measure at **one** peak
    point and scale analytically (voltage → nonlinear V–p curve; pulse length → PII;
    PRF → global I_spta multiplier).
- **Tiers:** Tier-0 anchor (1 pt/sequence) · Tier-1 validation cross per beam type ·
  bracket the worst case (max V × tightest focus for MI; max duty for I_spta).
- **RIS caveat:** the SFP uses **Random Interleaved Sampling** at fast timebases
  (`sample_mode=1`, `RIS_on` flag). RIS reconstructs one waveform from many repetitions —
  valid for repetitive imaging pulses, **invalid for single-shot ARF pushes**. Any rate
  > 100 MS/s on the 5112 is RIS. Use real-time single-shot for push/safety waveforms.
- **Record length / trigger:** capture window = `record_length / sample_rate` must cover
  the pulse; `reference_position` sets the pre/post-trigger split; set scope trigger
  **delay ≈ time-of-flight z/c** to place the window on the pulse arrival.
- **Saturation:** the preamp/scope path hits a hard ceiling at high drive. Use the
  unsaturated low/mid-voltage region and **extrapolate** linearly.

---

## 5. Calibration & MI computation

- **Scope input impedance:** a high-Z (1 MΩ) scope reads ~**2×** the 50 Ω voltage; the
  AH-2010 gain is referenced to a 50 Ω load. Apply **÷2 for 1 MΩ** (`LoadCorr`), ×1 for
  50 Ω. Verified internally consistent (50 Ω and 1 MΩ÷2 agree).
- **Parallel-circuit sensitivity:** `sens_total = G_amp · M_hyd · C_h/(C_h+C_amp)`
  ≈ **1.26 V/MPa** at 5 MHz.
- **Derating:** 0.3 dB/cm/MHz over the measured tip-to-face depth.
- **Acoustic working frequency f_awf:** IEC **−6 dB center of the fundamental**, from a
  **low-drive** capture (wideband centroid inflates with drive as harmonics build).
- **MI = pr.3 / √f_awf.**

---

## 6. L11-5 pulse-inversion MI validation (case study)

- **Sequence:** pulse inversion, transmit **1 cycle @ 4.5 MHz**. Received fundamental
  **f_awf ≈ 5.0 MHz** (−6 dB center; the L11-5 passband pulls the broadband 1-cycle
  pulse up). Verasonics likely uses 4.5 MHz in its MI calc.
- **Misaligned data** (elevation off): measured MI a constant **0.36×** Verasonics →
  diagnosed as **elevation misalignment** (constant multiplicative offset, shape
  preserved). Data archived in `MI L11-5\L11-5 misaligned\`.
- **Realigned data** (`MI L11-5\L11-5 realigned\`, z = 11 mm):
  - Shape linear and matches Verasonics; **50 Ω and 1 MΩ÷2 agree** (impedance ✓).
  - **Strong saturation** > ~35 V: 50 Ω ceiling ≈ 2.18 V, 1 MΩ ≈ 4.4 V (same TX voltage,
    2× apart → preamp-output ceiling). Use ≤ 30 V linear region + extrapolate.
  - Extrapolated MI ≈ **1.23×** Verasonics (at f_awf = 5.0), roughly constant.
  - Note files: filenames with `pos` = **normal** pulse, `neg` = **inverted** (notes in
    the `pos` files were corrected from "inverted" → "normal").
- **Residual ~1.23–1.3× — what it is NOT** (all checked):
  - _Frequency convention:_ using 4.5 MHz → ratio **1.32** (worse; MI ∝ 1/√f).
  - _Spatial averaging (400 µm tip):_ CF ≈ 1.03–1.07 and **raises** MI (wrong way).
  - _Broadband magnitude+phase deconvolution:_ +~2% (up to +11% at max drive), **raises**
    MI (wrong way).
  - _Reflections:_ echo at 2·11 mm/c ≈ 14.7 µs, outside the 5 µs capture window.
  - _Bubbles:_ attenuate (would read low) and are erratic; error is a steady over-read.
  - _Temperature:_ MI doesn't use ρc; sensitivity temp-coeff is sub-dB.
  - _Derating depth:_ 11 mm is a confirmed true tip-to-face distance.
- **Conclusion:** the residual **~1.3× (≈2.4 dB)** is within combined hydrophone-vs-system
  MI uncertainty. Remaining candidates: **absolute calibration accuracy (~±1–2 dB)** and/or
  **Verasonics' model reading low**. Caveat for write-up: measured **>** Verasonics, i.e.
  the non-conservative direction.
- **Old uncorrected 1 MΩ method = 2× overestimate** (missing ÷2 load correction).
- Figures: `MI L11-5\L11-5 realigned\MI_vs_TXvoltage_L11_5_realigned.png` (+ `_broadband`).

---

## 7. Beam-geometry reference & S5-1 example

For a focused aperture of width `D` at focal depth `z`, frequency `f`, sound speed `c`:

```
F# = z / D
lambda = c / f
lateral FWHM (-6 dB)  ≈ F# · lambda
depth of focus (-6 dB) ≈ 7.1 · F#² · lambda     (approximate)
derating factor        = 10^(0.3 · z[cm] · f[MHz] / 20)
```

**S5-1, full aperture (80 el × 0.254 mm), focus 79 mm, f = 1.95 MHz, c = 1480 m/s:**

| Quantity | Value |
|---|---|
| aperture D | 20.32 mm |
| F# | 3.89 |
| λ | 0.759 mm |
| lateral FWHM | ≈ 2.95 mm |
| **depth of focus (−6 dB)** | **≈ 82 mm** |
| derate @ focus (79 mm) | 4.6 dB (×1.70) |

**Axial safety sweep (per §4 recommendation):**
- **Span ≈ 40–118 mm** (0.5–1.5× focus) — covers the ~82 mm −6 dB DOF **and extends
  proximal**, because derating grows with depth and pulls the derated-MI/I_spta max
  _proximal_ of the 79 mm in-water focus.
- **Spacing:** two-stage — coarse **~5 mm** to map (~16 pts), then **~2 mm** refine within
  ±8–10 mm of the observed max. (Single-pass alternative: uniform ~3 mm.)
- At each axial point, run a small **transverse grid search** (lateral × elevation),
  span ≈ ±3 mm (≈ ±1 FWHM) at ~0.6 mm steps, save only the local-peak waveform + its
  coordinates. One sweep serves both MI and I_spta (both derated maxima track together).

---

## 7a. S5-1 focused transmit — results (COMPLETE)

- **Peak location:** `foc_6` = **(−77.41, 50.14, −46.60) mm**, depth ≈ 73 mm — the
  **derated** max (proximal/shallow edge of the broad in-water plateau; derating pulls it
  ~9 mm proximal of the plateau centre). Axial sweep 15 pts × 3 mm confirmed a bracketed
  max. Data: `S5-1\S5-1 10V focused\` (+ `S5-1_axial_profile.png`).
- **f_awf = 2.10 MHz** (measured −6 dB centre; programmed transmit 1.95 MHz).
- **Voltage sweep** (`S5-1\S5-1 voltage sweep focused\`, 2–20 V, 50 Ω): linear region
  4–12 V; saturates > ~15 V (preamp ceiling ~2.2 V). Extrapolate from the linear region.
- **Scanned mode:** PRI 270 µs × **71 focus lines** → frame period 19.17 ms →
  **local PRF at any point = 52.2 Hz** (each spot hit once/frame). This makes I_spta small.
- **Safe-voltage limits (linear extrapolation, depth 73 mm):**
  | Limit | Crosses at |
  |---|---|
  | MI = 1.9 | ~50 V |
  | I_sppa.3 = 190 W/cm² | ~48 V (binding) |
  | I_spta.3 = 720 mW/cm² | ~450 V (non-binding — scanned) |
- **Max safe TX ≈ 48 V (I_sppa.3-limited)**; recommend a working cap of **~45 V** given the
  ~4× extrapolation and depth/f_awf uncertainty. Plot: `S5-1_focused_safe_voltage.png`.
- Pulse duration ≈ 0.8 µs; only `pos` (normal) polarity measured.

## 7b. S5-1 push + full safety table (2026-08-14)

Depth for derating now taken from **`Home` = transducer centre** (= 37.5 mm at the push peak,
55 mm at the imaging peak). Push duty = **20 Hz for 1.2 s, then ≥30 s off → 0.77 Hz effective**.

| Configuration | TX (V) | MI | I_sppa.3 (W/cm²) | I_spta.3 (mW/cm²) |
|---|---|---|---|---|
| Focused imaging | 30 | 1.33 | 104 | 4.5 |
| Push 79 el, 1900 cyc, max V | 20.0 | 1.74 | 190 (bind) | 125 |
| Push 61 el, 1900 cyc, max V | 24.3 | 1.71 | 190 (bind) | 124 |
| Push 61 el, 1500 cyc, max V | 24.3 | 1.71 | 190 (bind) | 98 |
| Push 41 el, 1500 cyc, 30 V (est) | 30 | 1.52 | 144 | 75 |
| **FDA limit** | — | **1.9** | **190** | **720** |

- **Pushes are I_sppa.3-limited** (MI ~1.7, I_spta ~100–125 far under 720 thanks to the burst).
- **Pulse length scales I_spta only** — `ConfirmPulseLength.m`: MI ratio 1.01, I_sppa 0.98,
  I_spta 0.77 ≈ 1500/1900.
- **Focal gain ~1/F# (pressure ~ N^0.6–0.8):** fewer elements → higher allowed V, weaker/wider push.
- **41 el is estimated** (not measured; MI ~1.5–1.8, I_sppa ~144–199).
- Reproducible code: `SafetyTable.m`, `GenerateSafetyFigures.m`, `ConfirmPulseLength.m` (in `S5-1\`).

**TODO (measure to validate extrapolations):** (1) **25 V, 61 el, 1900 cyc** (recommended);
(2) **30 V, 41 el, 1900 cyc** (alternative — confirm MI < 1.9 and I_sppa). Also: transient TI for the
1.2 s burst; L12-3 safety table.

## 7c. 79-element no-preamp repeat variance (2026-08-18) — `RepeatVariance.m`

Each `.hws` in `S5-1\...\79 elements repeats\` holds **three** waveform records
(`/wfm_group0..2`, three consecutive push firings; the default `readHWS` reads only group 0).
Per-firing Vpp / Vneg vs the operator's multi-second visual observation:

| TX (V) | Vpp mean / hyp (mV) | Vneg mean / hyp (mV) | note |
|---|---|---|---|
| 30 | 103 / 102 | 51 / 52 | tight (103,103,103) |
| 40 | 130 / 129 | 65 / 66 | tight (130,129,130) |
| 45 | 142 / 142 | 72 / 74 | tight |
| 50 | 166 / 149 | 79 / 81 | **spread** (Vpp 155–183) |

- **Confirmed:** the mean trend matches the hypothesis within ~1–2 mV through 45 V, and is
  **monotonic** — no collapse. The earlier single-shot "dips" at 45–50 V were the supply-sag
  variance, caught at a low-transmit moment, not a real saturation drop.
- Only 50 V shows genuine firing-to-firing spread (worst supply sag at the highest 79-el load).
- Confirmed by the phantom Readme's Verasonics notice: **est. HV load 0.6 A > 0.5 A supply → actual
  transmit may drop below the selected value.**

## 7d. Frequency-domain vs fixed-f0 sensitivity (2026-08-18) — `FreqDomainSensitivity.m`

The default chain converts voltage→pressure with the loaded sensitivity at the fundamental only
(`p(t)=y/S(f0)`). As TX voltage rises the push develops stronger harmonics (2nd/1st ratio grows
0.25→0.44 over 20→50 V). Deconvolving in the frequency domain — dividing each bin by the
**frequency-dependent** loaded sensitivity `S(f)=M_oc(f)·C_h(f)/(C_h(f)+C_load)` across 1–12 MHz —
raises the indices modestly and in the **conservative** direction:

| TX (V) | H2/H1 | MI fix → fd | I_sppa fix → fd | ΔPr |
|---|---|---|---|---|
| 30 | 0.35 | 1.42 → 1.52 | 91 → 93 | +7% |
| 40 | 0.39 | 1.84 → 1.98 | 129 → 130 | +8% |
| 50 | 0.44 | 2.32 → 2.34 | 161 → 160 | +1% |

- Effect is **up to ~8 % on pr/MI**, ~1 % on I_sppa (the intensity integral is dominated by the
  fundamental). Direction is non-conservative→conservative, so ignoring it slightly **under**-reports.
- **Phase CAVEAT:** a full complex deconvolution needs the hydrophone **phase** vs frequency. The
  Onda HGL-0400 calibration is **magnitude-only**, so `S(f)` here is zero-phase and only the harmonic
  *amplitudes* are corrected (the hydrophone's own phase distortion is not). Treat this as a
  magnitude-only refinement until a phase-calibrated sensitivity is obtained. The AH-2010 preamp
  *does* have a calibrated phase (`AmplifierTable.PHASE_DEG`), used by `voltageToPressureBroadband.m`.
- **Dataset CAVEAT:** the absolute MI/Isppa above are computed on the **`79 elements repeats`**
  (typical, lower, monotonic) captures, so they are a fair basis for the *relative* fix-vs-fd
  comparison but **NOT for the safety limit**. The limit uses the worst-case peak-optimised
  `_opt_79el.hws` captures (Vneg ~1.6× higher) — see §7e.

## 7e. Max-voltage table CONFIRMED by direct no-preamp measurement (2026-08-18) — `CorrectedMaxVoltage.m`

Independent cross-check of the `SafetyTable` max voltages, which were obtained by extrapolating a
linear/quadratic fit of the low-voltage 2–8 V **with-preamp** sweep. Reading the **direct no-preamp**
peak-optimised push captures (15–50 V, `_opt` files, 41/61/79 el) and evaluating MI.3/Isppa.3/Ispta.3
at each measured voltage reproduces the table almost exactly:

| Config | MI=1.9 | Isppa=190 | Ispta=720 (1500/1900) | **max V** | table |
|---|---|---|---|---|---|
| 41 el | 39.3 V | **32.0 V** | >50 V | **32 V** | 32 |
| 61 el | 27.9 V | **23.1 V** | >50 V | **23 V** | 23 |
| 79 el | 23.1 V | **20.1 V** | >50 V | **20 V** | 20 |

- **The table is accurate.** All configs are **Isppa.3-limited**; 1500 and 1900 cyc share the same
  max V (Ispta never binds). The low-V extrapolation is validated because the limits are crossed at
  20–32 V, **below** where nonlinear saturation / supply sag bite (79el Vneg peaks ~40 V then dips at
  50 V — above the limits).
- **Correction to an earlier note:** a transient claim that the table "over-read ~4× / true max ~48 V"
  was **wrong** — it used the lower-pressure `repeats` (typical) set, not the worst-case `_opt`
  captures. For a safety limit we use the worst case; the two datasets differ ~1.6× (peak vs typical)
  at the same position/depth (supply-sag variance), not a calibration error.

**Comparison to the LaTeX reference chain** (`S5-1\Ultrasound safety indices.txt`): our chain matches
it step-for-step (remove offset → /G → FFT → P(f)=V(f)/S(f) → IFFT → pr→MI; I=p²/ρc → PII →
I_sppa=PII/τ_p; I_spta=I_sppa·τ_p·PRF). The **only** substantive difference was fixed-f0 vs
frequency-domain `S(f)`; `FreqDomainSensitivity.m` now implements the frequency-domain version and
quantifies the gap (above). Remaining LaTeX items already handled elsewhere: electrical gain (÷G),
finite-aperture spatial averaging (`spatialAvgFactor.m`), 0.3 dB/cm/MHz derating, ρ/c. Outstanding
LaTeX item still open: hydrophone **phase** calibration (see caveat).

## 8. File / data locations

- Tools: `D:\Luuk van Knippenberg\Claude\MI estimation\` (+ `MI L11-5\`)
- L11-5 data: `MI L11-5\L11-5 realigned\`, `MI L11-5\L11-5 misaligned\`
- MI pipeline + calibration: `D:\Luuk van Knippenberg\Github\SWI\Mechanical index\`
- New (2026-08-18): `hydrophone_analysis\RepeatVariance.m`, `FreqDomainSensitivity.m`
- Phantom SW sweep (tasks 2/3, separate pipeline): `shearWaveProcessing\` +
  `docs\phantom_voltage_sweep.md`. Data: `...\2026_08_17 Phantom sweep elements cycles TXvoltage\`.
  **Long-path gotcha:** the OneDrive folder prefix + `AcquisitionParametersAndECG.mat` exceeds
  Windows MAX_PATH (260) so the Python (Win32) pipeline can't `stat` the `.mat` (Git Bash can, via
  POSIX). Work through a short **junction** (`mklink /J D:\swp_ph "<sweep folder>"`).
- MATLAB: R2025b (`-batch`); note script names cannot start with `_` (invalid identifier).
