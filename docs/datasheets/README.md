# datasheets

Reference documents for every instrument and standard in the measurement. All third-party
material, reproduced here so the repository is self-contained; **copyright remains with the
issuing organisation** (Onda Corporation, the US FDA, Verasonics).

## Hydrophone chain

Onda publishes these on its website; they are mirrored here so the
repository stays self-contained. Copyright remains with Onda Corporation.


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
| `FDA-Guidance-Diagnostic-Ultrasound-Systems-Transducers.pdf` | the FDA guidance that sets the limits used throughout: MI ≤ 1.9 and I<sub>sppa.3</sub> ≤ 190 W/cm², with an I<sub>spta.3</sub> ceiling that depends on the application (720 mW/cm² peripheral vessel, 430 cardiac, 94 fetal & other, 17 ophthalmic) and 0.3 dB/cm/MHz derating. US Government work — public domain, hence reproduced here in full. |

## Verasonics — not included

Verasonics ships an MI/intensity table for each of its stock scripts, as part of the Vantage
installation. **That is licensed customer documentation and is not redistributed here** — but you
already have it if you have a Vantage. It sits next to the script it describes, e.g.

```
Example_Scripts\Biomedical\Vantage 128 and 256\UTA-260-S and 260-D\
    L11-5v\High Image Quality\MIData_L11-5v-HIQScripts.pdf
```

with the same layout for the other probes (`C5-2v`, `P4-2v`, …). The `SetUp*.m` example scripts
this repository refers to are in those same folders.

What *is* here is our own reading of them: the thirteen MI values that `RunDemo` and `PlotMI_L11_5`
validate against are tabulated in the code, as measurements taken from the console for the L11-5v
`WideBeamHISC` sequence. That is data we recorded, not their document.

## Digitizer

| File | Why you would open it |
|---|---|
| `NI PCI-5112 specs.pdf` | the digitizer: **8-bit**, 100 MS/s real-time, software-selectable 50 Ω / 1 MΩ input. Three properties drive the protocol — 8 bits is why the vertical range must be set tightly, 100 MS/s is the ceiling above which the card switches to RIS (invalid for a single-shot push), and the switchable impedance is what makes the 50 Ω / 1 MΩ check a keystroke rather than a cable change. |

The Scope Soft Front Panel documentation is not bundled; it comes with the NI-SCOPE driver
installation. What matters for the analysis is in `readHWS`: `.hws` is HDF5, and samples are stored
as raw ADC codes that must be scaled by the polynomial in the file.
