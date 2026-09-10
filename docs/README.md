# docs

| Path | What |
|---|---|
| `report/` | **the protocol report** — how to measure the acoustic output of a custom sequence, written for someone who has to repeat it. Read [`report/main.pdf`](report/main.pdf); the source and figures are alongside it, and [report/README.md](report/README.md) says how to rebuild. |
| `HydrophoneSafety_Notes.md` | the working notes: setup, the `.hws` format, what was tried, what went wrong, what the findings were. Less polished than the report and more candid — read it when the report says *what* and you want *why*. |
| `UltrasoundSafetyIndices.pdf` | the processing chain written out as mathematics: V<sub>DAQ</sub>(t) → V<sub>hyd</sub>(f) → P(f) → p(t) → {MI, I<sub>sppa</sub>, I<sub>spta</sub>} — the general frequency-domain form of what `SafetyIndices.m` evaluates at a single frequency and `FreqDomainSensitivity.m` evaluates in full. |
| `UltrasoundSafetyIndices.tex` | its source. Standalone: `tectonic -X compile UltrasoundSafetyIndices.tex --outdir .` |
| `datasheets/` | every instrument and standard involved — see [datasheets/README.md](datasheets/README.md). |

The calibration certificates for our two specific units are not here; they are with the
calibration data they describe, in [`../calibration/certificates/`](../calibration/certificates/).
