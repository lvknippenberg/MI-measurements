# calibration

The calibration of the hydrophone chain. Everything the code needs to turn volts into pascals.

| File | What |
|---|---|
| `HGL0400-1037_hydrophone_opencircuit_20240610.txt` | Onda HGL-0400 SN 1037: **open-circuit** end-of-cable sensitivity M<sub>c</sub>(f) and hydrophone capacitance C<sub>H</sub>(f), 1–20 MHz |
| `AH-2010-025-1109_amplifier_20240610.txt` | Onda AH-2010-025 SN 1109: gain G(f), phase, input capacitance C<sub>A</sub>(f), 50 Ω output |
| `HydrophoneTables.mat` | the two files above, digested into `HydrophoneTable` (FREQ_MHz, SENS_DB, CAP_PF) and `AmplifierTable` (FREQ_MHZ, GAIN_DB, PHASE_DEG, CAP_PF) — this is what the code loads, via `miCalibration` |
| `certificates/HGL0400-1037_hydrophone_certificate_20240610.pdf` | the signed Onda certificate as issued |
| `certificates/HGL0400-1037_hydrophone_certificate_20240610.txt` | its raw text export, as delivered |
| `certificates/AH-2010-025-1109_amplifier_certificate_20240610.pdf` | ditto, pre-amplifier |
| `certificates/AH-2010-025-1109_amplifier_certificate_20240610.txt` | ditto, pre-amplifier |

Both instruments were calibrated by Onda on **2024-06-10**.

## The one thing to get right

The hydrophone certificate says `ELECTRICAL_LOAD OpenCircuit`. An open-circuit sensitivity is
**not** what the scope sees. Converting it to the loaded sensitivity requires the capacitive
divider between the hydrophone and whatever loads it (Onda `HydroCalMethod` Eq. 2a — the PDF is
in [`../docs/datasheets/`](../docs/datasheets/)):

```
with the preamp:      sens = G(f) · M_c(f) · C_H / (C_H + C_A)        C_A ≈ 6.1 pF
without the preamp:   sens =        M_c(f) · C_H / (C_H + C_load)     C_load ≈ 123 pF
```

Dropping the divider under-reads the pressure by 1.5×.

`C_load` — the cable, the adapters and the scope input together — is not a specified quantity.
It is roughly ten times the hydrophone's own capacitance, so it dominates a no-preamp reading
entirely, and it has to be **measured**: see `CalibrateLoadCapacitance`, which recovers it from
paired with/without-preamp captures of the same acoustic point. Re-measure it whenever anything
in that path changes — a different cable, an added BNC tee, a different scope input.

## Regenerating `HydrophoneTables.mat`

The `.mat` is a convenience, not the source of truth: the two `.txt` certificates are. If you
recalibrate, replace the `.txt` files and rebuild the tables from their `FREQ`/`SENS`/`CAP`/`GAIN`
columns, keeping the variable and column names above — `miCalibration`, `SystemSensitivity` and
the analysis scripts all read them by name.
