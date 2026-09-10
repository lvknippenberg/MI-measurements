# demo

Two demonstrations. Start with the small one.

## `DemoSingleCapture` — one waveform in, the safety numbers out

```matlab
cd MI-measurements
misetup
DemoSingleCapture                      % a bundled example capture
DemoSingleCapture(myFile)              % any .hws capture
S = DemoSingleCapture(myFile, PRF_Hz=52.2)
```

The whole chain at its smallest: read a capture, work out how it was acquired from the note stored
inside it, convert volts to pressure, print MI / I<sub>sppa.3</sub> / I<sub>spta.3</sub> against the
limits, and write one figure showing where each number came from — the waveform with p<sub>r</sub>
and p<sub>+</sub> marked and the 10–90 % window shaded, a zoom on the rarefactional half-cycle, the
spectrum with the −6 dB band and f<sub>awf</sub>, and the cumulative intensity integral that fixes
the pulse duration.

It configures itself from the capture note via `CaptureConditions`: `Home`/`Pos` give the
depth, `Scope_impedance` picks the load correction, and `Pre-amp = No` switches to the
bare-hydrophone chain. **If the note is missing one of those it asks** — a menu for the
impedance, a number for the depth — and marks the value in its output as supplied rather than
recorded. Under `matlab -batch`, where nothing can be asked, it errors instead of guessing;
pass `Depth_cm`, `ScopeImpedance_Ohm` or `WithPreamp` to answer in advance.

Two further things it cannot know and takes as arguments:

- **`PRF_Hz`** — a property of the *sequence*, not of the capture. Without it I<sub>spta.3</sub> is
  reported as not evaluated rather than guessed. For a burst protocol pass the *effective* rate.
- **`LoadCap_pF`** — only for no-preamp captures; defaults to the 123.3 pF measured for this setup.
  Re-measure with `CalibrateLoadCapacitance` if the cabling changed.

## `RunDemo` — the full demonstration

```matlab
RunDemo
```

About a minute. Figures land in `results/demo/`. Neither demonstration needs anything outside
this repository, and neither suppresses warnings.

`RunDemo` walks the whole chain — calibration → one capture → a voltage sweep → validation
against Verasonics — and then checks every number it produced against
`reference/expected_results.csv`. The last line tells you whether the install is sound:

```
  ALL 22 QUANTITIES REPRODUCE THE REFERENCE.
```

## The reference file

`reference/expected_results.csv` is 22 numbers spanning all four parts: the calibration
constants, the indices of a single capture, the maximum allowed transmit voltage for each of the
three apertures, and the L11-5 validation ratio. They are compared at a relative tolerance of
1e-9 — the chain is fully deterministic, so anything looser would hide a real change.

Three of those numbers, the `sweep.el*.maxV` values, match `SafetyTableAll`'s output **to the
last digit**, which is the point: the compact `SafetyIndices` API used by the demo and the full
session-by-session analysis are the same computation, not two approximations of it.

Regenerate the reference only after a deliberate, understood change to the method:

```matlab
RunDemo('reference')
```

If a mismatch appears that you did *not* intend, find out why before regenerating — that check is
the only thing standing between a silent method change and a wrong safety limit.
