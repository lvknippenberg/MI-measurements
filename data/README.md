# data

Real hydrophone captures, bundled so every analysis in this repository runs out of the box.
193 `.hws` files, about 35 MB. Resolve paths with `miData(...)`, never by typing them.

`.hws` is NI **Hierarchical Waveform Storage** — HDF5 underneath, readable with `h5read` and no
NI drivers. `readHWS` handles the details: samples are stored as **raw 8-bit ADC codes** and must
go through the axis `scale_coef` polynomial to become volts.

## The note inside each capture

Every file carries the acquisition context as an HDF5 attribute, which is what makes the set
self-describing:

```
Probe = S5-1
Home = (55.8,105.84,-17.3) mm      <- transducer centre
Pos  = (46.1,70.84,-17.3) mm       <- hydrophone tip
TX_voltage = 50V
Pulse = normal
Scope_impedance = 1e6
Pre-amp = No
PushCycles = 1900
PushElements = 79
```

**`Home`, `Pos` and `Scope_impedance` are the three fields the analysis cannot do without**:
the first two give the derating depth `|Pos − Home|`, the third selects the load correction.
`CaptureConditions` prompts for whichever is absent rather than assuming one. The L11-5 set
logged `Pos` only, already relative to the probe — pass `Home_mm = [0 0 0]` for those. The
remaining fields are provenance and never enter the arithmetic.

## What is here

### `S5-1/2026-08-13_sessionA_preamp_50ohm/`

The shear-wave sequence measured **with** the AH-2010 pre-amplifier on a 50 Ω-terminated scope.
The preamp clips at ~2.2 V, so this session is confined to low transmit voltages and its limits
come from extrapolation.

| Folder | Files |
|---|---|
| `focused_imaging_sweep/` | focused imaging, 2–20 V at the peak |
| `focused_imaging_axial/` | focused imaging, 10 V, 15 axial stations |
| `push_79el_sweep/` | 79-element push, 2–12 V, plus one 1500-cycle capture at 5 V |
| `push_79el_axial/` | 79-element push, 5 V, 16 axial stations |
| `push_61el_sweep/` | 61-element push, 3–8 V |
| `push_61el_axial/` | 61-element push, 5 V, 13 axial stations |

### `S5-1/2026-08-17_sessionB_preamp_vs_nopreamp/`

The same push measured **without** the pre-amplifier, straight into the 1 MΩ input, so the whole
15–50 V range is measured directly with no clipping. Also holds the with/without pairs at 2–12 V
that overlap session A — those pairs are what `CalibrateLoadCapacitance` uses.

- `S5-1_<V>V_Preamp_peak.hws` — with preamp, at the peak, 2–20 V
- `S5-1_<V>V_NoPreamp_opt.hws` — without preamp, peak-optimised, 2–50 V (61-element base)
- `S5-1_<V>V_NoPreamp_opt_41el.hws`, `..._79el.hws` — the other two apertures
- `S5-1_4V_Preamp_<k>.hws` — a 10-station axial sweep at 4 V
- `79el_repeats/` — three consecutive firings stored per file, for the firing-to-firing spread

Two folders also carry a `Notes.txt` — the operator's acquisition note for that session, including the `Home` coordinate the analysis uses — and `push_79el_axial/` carries `S5-1_5V_push_coordinates.csv`, the stage position and peak-negative voltage recorded at each station of that sweep.

### `L11-5/2026-08-13_realigned/`

The L11-5v pulse-inversion sequence, 2–50 V, both polarities, on both 50 Ω and 1 MΩ inputs.
This is the **validation** set: Verasonics reports an MI for this sequence, so the whole chain
can be checked against an independent number. ("Realigned" — an earlier set taken before the
elevation alignment was corrected read low and is not included.)

## Two things to know before trusting a number from these files

**Peak-optimised vs. repeat captures.** Two 79-element datasets exist at the same nominal
position and they differ by ~1.6×. The three firings *inside* each repeat capture agree to 1–2 %,
so this is not supply sag — it is a between-capture positioning difference, and the repeats sample
slightly off-peak. Session A, an independent measurement on a different day with a different
chain, agrees with the **peak-optimised** captures to ~8 %. Those set the limits; the repeats must
not be used for one.

**Saturation.** With the preamp, captures above roughly 10 V (push) or 25 V (L11-5) are clipped
by the AH-2010, not by the field. The analyses fit the unsaturated region and extrapolate; a raw
reading from a clipped capture is meaningless and cannot be rescued by rescaling.

## Using a different data set

`miData` honours the `MI_DATA_ROOT` environment variable. Point it at a share holding the full,
unabridged measurement set and every script follows, with no edits:

```matlab
setenv('MI_DATA_ROOT', '\\server\share\hydrophone')
```
