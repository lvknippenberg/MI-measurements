# MI-measurements

Measuring the **acoustic output** — mechanical index MI, and the intensities
I<sub>sppa.3</sub> and I<sub>spta.3</sub> — of a custom Verasonics transmit sequence with a
calibrated hydrophone, and deciding from that what transmit voltage the sequence may use.

Everything needed is in this repository: the measurement protocol, the calibration of the
hydrophone chain, the datasheets for every instrument, stand-alone MATLAB code, real example
captures, and a demonstration that runs end to end and checks itself against stored reference
numbers.

The measurements documented here were made for a cardiac **shear-wave elastography** sequence
(acoustic-radiation-force push + tracking) on a Verasonics Vantage with an S5-1 phased array,
with an L11-5v pulse-inversion sequence used as an independent cross-check against the MI the
console reports.

---

## Quick start

```matlab
cd MI-measurements
misetup              % puts code/, code/analysis/, code/acquisition/, demo/ on the path
DemoSingleCapture    % one capture -> the safety numbers and an annotated figure
RunDemo              % the full demonstration, ~1 min, writes figures to results/demo/
```

[`DemoSingleCapture`](demo/DemoSingleCapture.m) is the chain at its smallest and the place to
start: give it any `.hws` file and it prints MI, I<sub>sppa.3</sub> and I<sub>spta.3</sub> against
the limits and draws the waveform with every quantity it read off it. It works out the measuring
chain from the note stored inside the capture.

`RunDemo` needs nothing but MATLAB (developed on R2025b) and this repository — no Verasonics,
no NI drivers, no absolute paths. It finishes with a **regression check**: every number it
computed is compared against `demo/reference/expected_results.csv`, so a broken install or a
changed method shows up immediately rather than silently.

What the five parts of the demonstration do:

| Part | What it shows |
|---|---|
| 1 Calibration | the two Onda certificates → a sensitivity in V/MPa, for the chain with and without the pre-amplifier (including the cable+scope load capacitance, which is *measured*, not looked up) |
| 2 One capture | a single `.hws` file → p<sub>r</sub>, f<sub>awf</sub>, pulse duration, PII, MI, I<sub>sppa.3</sub>, I<sub>spta.3</sub>, with a figure showing where each number comes from |
| 3 Safety limit | a transmit-voltage sweep for 41/61/79-element pushes → the maximum voltage at which the sequence still meets the FDA limits, and which of the three binds |
| 4 Validation | the same chain applied to the L11-5 sequence, compared against the MI Verasonics reports for it |
| 5 Regression | all of the above vs. the stored reference |

---

## What is in here

```
misetup.m              put the repository on the MATLAB path
code/                  the measurement chain, stand-alone
  analysis/            the full analyses behind the report's figures and tables
  acquisition/         scan planning and the Verasonics safety-evaluation sequences
  legacy/              the pre-2026-08 code, kept for provenance -- do not use
calibration/           the Onda certificates and the digested tables
data/                  real example captures (193 .hws files, ~35 MB)
demo/                  DemoSingleCapture.m (one capture), RunDemo.m (everything) + reference
docs/                  the protocol report (built PDF + source), working notes, all datasheets
results/               anything the code writes (git-ignored)
```

### The core API (`code/`)

Five functions are the whole measurement chain. Everything else is built on them.

| Function | What it does |
|---|---|
| [`readHWS`](code/readHWS.m) | read an NI **Hierarchical Waveform Storage** capture (`.hws` is HDF5) → waveform, full scope configuration, and the note written at acquisition time. Warns on RIS (equivalent-time) sampling and on a record shorter than the pulse. |
| [`AcousticWorkingFrequency`](code/AcousticWorkingFrequency.m) | f<sub>awf</sub> as IEC 62127 defines it: the centre of the −6 dB band of the fundamental. Measure it on a *low-drive* capture — harmonics bias it upward. |
| [`SystemSensitivity`](code/SystemSensitivity.m) | calibration → V/MPa, for either chain, including the scope-impedance correction |
| [`CalibrateLoadCapacitance`](code/CalibrateLoadCapacitance.m) | the cable+scope capacitance a bare hydrophone sees, from paired with/without-preamp captures |
| [`SafetyIndices`](code/SafetyIndices.m) | waveform → p<sub>r</sub>, MI, PD, PII, I<sub>sppa.3</sub>, I<sub>spta.3</sub> |

[`CalculateSafety`](code/CalculateSafety.m) is the older converter that goes the other way —
from a single measured V<sub>min</sub> plus a Verasonics `TW` structure whose pulse shape it
*simulates*. It is kept because the Verasonics-side workflow uses it, but it needs
`computeTWWaveform` from Verasonics on the path (or a `TW` that already carries `Wvfm2Wy`).
**Prefer `SafetyIndices`**: it uses the pulse that was actually measured, so nothing depends on
a transmit model.

### The analyses (`code/analysis/`)

| Script | What it produces |
|---|---|
| `SafetyTableAll` | the authoritative table: every S5-1 session, one set of conventions, maximum voltage per aperture, plus the cross-validation between the two independent methods |
| `SafetyTable` | the earlier session-A-only table (extrapolated); superseded, kept for comparison |
| `CorrectedMaxVoltage` | session-B-only cross-check of the same limits |
| `PreampComparison` | with vs. without the pre-amplifier: do they agree, and what is C<sub>load</sub> |
| `MakeReportFigures` | the four figures the report builds from the captures, plus the tank-reflection diagnostic |
| `GenerateSafetyFigures` | MI vs. axial distance and MI vs. transmit voltage, per configuration |
| `ConfirmPulseLength` | direct data confirmation that pulse length scales I<sub>spta</sub> only |
| `PlotMI_L11_5` | the L11-5 vs. Verasonics comparison, with saturation-aware extrapolation |
| `RepeatVariance` | firing-to-firing spread inside one capture (transmit supply sag) |
| `FreqDomainSensitivity` | fixed-f<sub>0</sub> vs. full frequency-domain deconvolution |
| `PlotPushVoltage` | peak-to-peak and peak-negative voltage vs. drive, with the saturation zone |

They write to `results/`, never back into `data/`.

### The protocol report (`docs/report/`)

**[`docs/report/main.pdf`](docs/report/main.pdf)** is the document to hand a colleague who has to
repeat this on another sequence: the setup, the alignment procedure, the calibration mathematics,
the pitfalls, and the worked results. It is checked in, so no TeX toolchain is needed to read it;
the source is `main.tex` alongside it. To rebuild after editing the source, see
[docs/report/README.md](docs/report/README.md) — Tectonic locally (a single self-contained
executable) or Overleaf.

---

## The measurement setup

| | |
|---|---|
| Hydrophone | Onda **HGL-0400**, SN 1037, 400 µm aperture, calibrated 2024-06-10 |
| Pre-amplifier | Onda **AH-2010-025**, SN 1109, 20 dB, 6.1 pF input, **50 Ω output** |
| Digitizer | NI **PCI-5112**, 8-bit, 100 MS/s real-time, software-selectable 50 Ω / 1 MΩ input |
| Software | NI Scope Soft Front Panel 3.8, writing `.hws` (HDF5) |
| Positioner | 3-axis stage; depth is `|Pos − Home|` from the note stored in each capture |
| Ultrasound | Verasonics Vantage, S5-1 phased array (and L11-5v for the cross-check) |

Datasheets and manuals for all of it are in [`docs/datasheets/`](docs/datasheets/), the
calibration certificates in [`calibration/certificates/`](calibration/certificates/), and the
FDA guidance that sets the limits is
[`docs/datasheets/FDA-Guidance-Diagnostic-Ultrasound-Systems-Transducers.pdf`](docs/datasheets/).

---

## Two corrections that matter

Both were found in August 2026 and both are already applied throughout this repository. They are
worth stating plainly, because the first one made every earlier number wrong by a factor of two.

**1 — Scope input impedance (MI ×2, intensities ×4).**
The AH-2010 has a **50 Ω output impedance**, and its 20 dB gain — hence the Onda combined system
sensitivity — is referenced to a **50 Ω load**. Driving a high-impedance (1 MΩ) scope, it delivers
about twice the voltage it delivers into 50 Ω. Every pressure read on a 1 MΩ input was therefore
~2× high, and every intensity ~4× high. The fingerprint was visible in the data all along: the
preamp clips at 4 V<sub>pp</sub> (2 V peak) into 50 Ω, yet the observed "saturation" ceiling was
~4.8 V peak.

*Either terminate the scope into 50 Ω (best), or pass the true input impedance and let
`SystemSensitivity` absorb the factor.*

**2 — The capacitive divider is not optional.**
The Onda certificate is an **open-circuit** sensitivity (`ELECTRICAL_LOAD OpenCircuit`), so
converting it to a loaded sensitivity *requires* the divider term C<sub>H</sub>/(C<sub>H</sub>+C<sub>A</sub>)
≈ 0.667 (Onda `HydroCalMethod` Eq. 2a). Omitting it under-reads the pressure by 1.5×. This is the
`"ParallelCircuit"` method, and it is the default everywhere here; the older `"Simple"` method is
wrong.

Without the preamp the same divider applies, but against the **cable+scope** capacitance instead —
about 123 pF against the hydrophone's 12 pF, so it dominates the reading completely. That number
is on no datasheet. It must be measured, which is what `CalibrateLoadCapacitance` does, and it
must be re-measured whenever the cable, an adapter, or the scope input changes.

---

## Results for the S5-1 shear-wave sequence

From `SafetyTableAll` (and reproduced by `RunDemo`). Push duty is 24 pushes over a 1.2 s burst
followed by ≥30 s idle → an effective 0.77 Hz, which is what I<sub>spta.3</sub> uses.

| Configuration | max TX (V) | MI there | I<sub>sppa.3</sub> there | Limited by |
|---|---|---|---|---|
| Push, 41 elements | 32.0 | 1.48 | 190 | I<sub>sppa.3</sub> |
| Push, 61 elements | 23.1 | 1.47 | 190 | I<sub>sppa.3</sub> |
| Push, 79 elements | 20.1 | 1.62 | 190 | I<sub>sppa.3</sub> |
| *FDA limits (peripheral vessel)* | — | *1.9* | *190 W/cm²* | — |

Three things fall out of this:

- **The pushes are I<sub>sppa.3</sub>-limited**, never MI-limited. MI has headroom left at the
  voltage where I<sub>sppa.3</sub> binds.
- **I<sub>spta.3</sub> is nowhere near its limit** — because of the burst duty. Under a
  *continuous* 20 Hz assumption it would be ~26× higher and would bind first, at around 10 V. The
  burst is what buys the headroom, so the burst is part of the safety case. (The transient
  temperature rise across each 1.2 s burst is a separate question, to be checked via TI.)
- **Pulse length scales I<sub>spta</sub> only.** MI and I<sub>sppa</sub> are unchanged —
  confirmed directly from 1900- vs 1500-cycle captures (`ConfirmPulseLength`: MI ratio 1.01,
  I<sub>sppa</sub> 0.98, I<sub>spta</sub> 0.77 ≈ 1500/1900).

The two sessions reach these limits by **independent routes** — session A extrapolates a
low-voltage sweep taken *with* the preamp on a 50 Ω input, session B measures the full 15–50 V
range directly *without* the preamp on 1 MΩ — and they agree to 2–7 %. That agreement is the main
reason to trust the numbers.

Against Verasonics, on the L11-5 pulse-inversion sequence, the corrected chain gives
MI ≈ **1.23×** the console value, within the combined uncertainty of a hydrophone MI. The
uncorrected 1 MΩ reading gave ~2× — which is how correction 1 above was found.

---

## Open items

- The **41-element row rests on session B alone.** Session A has no 41-element data, so it has
  not been cross-validated the way 61 and 79 have. Measure it before relying on it.
- **Confirm the ×2 impedance factor directly** with a 1 MΩ-vs-50 Ω terminator measurement on an
  unsaturated signal. The evidence for it so far is indirect (the clip ceiling, and the agreement
  of the two sessions).
- **I<sub>spta.3</sub> and probe heating** need recomputing before any *continuous* high-PRF
  transmit — for imaging or passive elastography, that budget binds long before MI does.
- The report keeps the **ideal** measurement (a full 3-D raster with the Verasonics, stage and
  digitizer synchronised) in §2.1 alongside the reduced axial-line measurement actually
  performed. If that synchronisation ever gets built, §2.1 is the specification to build against.

---

## Provenance

The code and the report were developed in the
[`SWI`](https://github.com/lvknippenberg/SWI) sequence repository (`Mechanical index/`), with the
raw captures held outside it. This repository is the stand-alone consolidation: same code, same
calibration, now with the data alongside it and every path resolved relative to the repository
root. `SafetyTableAll` here reproduces the original output **bit-for-bit**.

Working notes with the fuller narrative are in
[`docs/HydrophoneSafety_Notes.md`](docs/HydrophoneSafety_Notes.md); the processing-chain
derivation is in [`docs/UltrasoundSafetyIndices.tex`](docs/UltrasoundSafetyIndices.tex).
