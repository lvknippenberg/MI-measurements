# code

MATLAB. Developed on R2025b, base MATLAB throughout the core; the analysis scripts additionally
use `hann` and `db2mag` from the Signal Processing Toolbox. Run `misetup` once per session — every
path resolves relative to the repository root, so nothing here needs editing after a clone.

## The chain (`code/`)

These five functions are the measurement. Read them in this order.

```matlab
[w, cfg] = readHWS(file);                                  % capture -> waveform + config + note
C        = CaptureConditions(cfg);                         % depth, chain, scope impedance
fawf     = AcousticWorkingFrequency(w(1).y, w(1).dt, [1.3 3.5]);
sens     = SystemSensitivity(fawf, Chain=C.chain, ScopeImpedance_Ohm=C.scopeImpedance_Ohm);
S        = SafetyIndices(w(1).y, w(1).dt, Sensitivity=sens, ...
                         Depth_cm=C.depth_cm, Fawf_MHz=fawf, PRF_Hz=0.769);
```

`DemoSingleCapture` is exactly this, with a figure and a printed summary.

| File | |
|---|---|
| `readHWS.m` | `.hws` (HDF5) → waveform, scope configuration, acquisition note. Warns on RIS sampling and on a record shorter than the pulse. |
| `quietHDF5TimestampWarning.m` | silences the one unavoidable HDF5 datatype warning (NI 128-bit timestamps) by identifier, around `h5info` only. Never suppress warnings globally instead — the other two from `readHWS` are real. |
| `AcousticWorkingFrequency.m` | f<sub>awf</sub>: centre of the −6 dB band of the fundamental. Also accepts a filename. |
| `SystemSensitivity.m` | V/MPa for either chain, with the scope-impedance correction folded in. |
| `CalibrateLoadCapacitance.m` | C<sub>load</sub> (cable + adapters + scope) from paired with/without-preamp captures. |
| `DepthFromNotes.m` | `\|Pos − Home\|` from the capture note. |
| `CaptureConditions.m` | resolves the acquisition conditions for one capture. **Home, Pos and Scope_impedance are required**; when the note lacks one it prompts for it interactively, errors under `-batch` rather than guessing, and labels every value with where it came from. `Pre-amp` is the one field with a default (absent = it was in the chain). |
| `SafetyIndices.m` | waveform → p<sub>r</sub>, MI, PD, PII, I<sub>sppa.3</sub>, I<sub>spta.3</sub>. |
| `voltageToPressureBroadband.m` | full frequency-domain deconvolution (hydrophone magnitude + amplifier phase), for harmonic-rich waveforms. |
| `spatialAvgFactor.m` | finite-aperture correction, jinc focal-plane model. First-order only — see its caveats. |
| `CalculateSafety.m` | the older Verasonics-side converter: one V<sub>min</sub> + a `TW` structure → the indices, **simulating** the pulse. Needs `computeTWWaveform` from Verasonics unless `TW.Wvfm2Wy` is already present. Prefer `SafetyIndices`. |
| `miRoot.m`, `miData.m`, `miCalibration.m` | path and calibration resolution. |

## `analysis/`

The full analyses. All of them write to `results/`, never back into `data/`.
`SafetyTableAll` is the authoritative one — it processes every S5-1 session on one set of
conventions and supersedes `SafetyTable` (session A only) and `CorrectedMaxVoltage` (session B
only), both of which are kept so the two independent routes can still be compared.

## `acquisition/`

Everything used at the water tank rather than afterwards. Captures themselves are taken by hand
in the NI Scope Soft Front Panel, one file per point, with the acquisition context typed into
the note field — that note is what `DepthFromNotes` later reads.

| File | |
|---|---|
| `scanPlan.m` | a non-uniform 3-line scan sized to the beam (λ, lateral FWHM, depth of focus), with per-point trigger delays. |
| `measurementPlan.m` | the per-beam-type protocol: coarse search → recentre on the observed peak → fine scan → 1-D voltage and pulse-length sweeps at that peak. |

No transmit sequences are included: they are Verasonics code, or derived from it, and cannot be
redistributed. All three were the same reduction of a stock example — steering removed, a single
transmit event repeated continuously so the scope can self-trigger (\S"Step 1" of the report):

| Sequence | Derived from |
|---|---|
| S5-1 shear-wave push and imaging (produced everything in `data/S5-1/`) | our `SetUp_SWI_Widebeam` |
| L11-5v pulse inversion (produced `data/L11-5/`) | `SetUpL11_5vWideBeamHISC.m` |
| L12-3v | `SetUpL12_3vWideBeamSC.m` |

The Verasonics originals are under
`Example_Scripts\Biomedical\Vantage 128 and 256\UTA-260-S and 260-D\<probe>\` in the Vantage
installation.

## `legacy/`

Two files, kept only so the corrections are traceable. **Do not use them**: they predate both the
scope-impedance correction and the `ParallelCircuit` default, and they contain hardcoded paths to
folders that no longer exist.

| File | |
|---|---|
| `Legacy_CalculateSafety.m` | what produced the numbers later found to be ~2× high in MI and ~4× in intensity |
| `Legacy_CalculateMI.m` | the driver script that called it |
