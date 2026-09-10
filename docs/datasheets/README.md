# datasheets

Reference documents for every instrument and standard in the measurement. All third-party
material, reproduced here so the repository is self-contained; **copyright remains with the
issuing organisation** (Onda Corporation, the US FDA, Verasonics).

## Hydrophone chain

| File | Why you would open it |
|---|---|
| `Onda_HGL_DataSheet.pdf` | HGL-series specifications: aperture, frequency range, nominal sensitivity, capacitance, maximum pressure |
| `Onda_AH-2010_DataSheet.pdf` | the pre-amplifier: 20 dB gain, 25 MHz bandwidth, 6.1 pF input capacitance, **50 Ω output**, and the ~4 V<sub>pp</sub> output ceiling that clips every high-drive push measurement |
| `Onda_HydroCalMethod.pdf` | **the important one.** Onda's calibration method, including Eq. 2a — the open-circuit-to-loaded conversion with the capacitive divider that must not be dropped. |
| `OndaCombineCal20151210.pdf` | the OndaCombineCal utility, which combines hydrophone and amplifier calibrations into a system sensitivity |
| `Onda_HydrophoneCare.pdf` | handling, cleaning, storage. A 400 µm needle tip is easy to destroy and expensive to replace. |
| `Onda_Hydrophone_Handbook_Link.pdf` | pointer to Onda's fuller hydrophone handbook |

## Standard

| File | |
|---|---|
| `FDA-Guidance-Diagnostic-Ultrasound-Systems-Transducers.pdf` | the FDA guidance that sets the limits used throughout: **Track 3** — MI ≤ 1.9, I<sub>sppa.3</sub> ≤ 190 W/cm², I<sub>spta.3</sub> ≤ 720 mW/cm², with 0.3 dB/cm/MHz derating. Also defines the application-specific I<sub>spta.3</sub> ceilings. |

## Verasonics

The MI/intensity data the console reports for its stock scripts, per probe. `Verasonics_MIData_L11-5v.pdf`
is the source of the reference MI curve that part 4 of `RunDemo` validates against.

| File | Probe |
|---|---|
| `Verasonics_MIData_L11-5v.pdf` | L11-5v |
| `Verasonics_MIData_C5-2v.pdf` | C5-2v |
| `Verasonics_MIData_P4-2v.pdf` | P4-2v |
| `Verasonics_AcousticSafety_PULSE_example.pdf` | worked example of the acoustic-safety calculation on a Vantage |

## Not included

The **NI PCI-5112** manual and the Scope Soft Front Panel documentation are not bundled — they
come with the NI-SCOPE driver installation. What matters for the analysis is captured in
`../HydrophoneSafety_Notes.md` and in `readHWS`: 8-bit, 100 MS/s real-time, software-selectable
50 Ω / 1 MΩ input, and `.hws` files that are HDF5 with samples stored as raw ADC codes.
