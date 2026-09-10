# The MI / acoustic-safety measurement protocol

`main.tex` is the document: how to measure **MI, I<sub>sppa.3</sub> and I<sub>spta.3</sub>** of a
custom Verasonics sequence with an Onda hydrophone and an NI PCI-5112. It is written as a protocol
for a colleague to follow, not as a log of what we did — our own numbers appear only in the
"Worked results" section.

| Path | What |
|---|---|
| `main.pdf` | **the built report** -- read this; it is checked in so the repository is useful without a TeX toolchain |
| `main.tex` | the source |
| `flowchart.tex` | the one-page overview figure (TikZ), `\input` at the head of the protocol section |
| `owis.tex` | Appendix A, operating the OWIS motion stage. Adapted from the standalone note *OWIS motion stage -- how to setup and use* (19 Mar 2024); its five figures were extracted from that PDF into `figures/owis_*`. **The OWISoft installation password in the original is deliberately not reproduced** -- do not paste it back in. |
| `nomenclature.tex` | the symbols-and-abbreviations list, `\input` after the table of contents. Unnumbered on purpose: making it a numbered section would shift every section number, and the code and READMEs cite those. |
| `figures/` | every figure it includes, including the two setup photographs |

The code that regenerates the figures lives with the rest of the code, in
[`../../code/analysis/`](../../code/analysis/) — it is not duplicated here.

## Building the PDF

**Locally, with Tectonic** (recommended -- one 20 MB executable, no admin rights, no TeX Live
install; it downloads the packages the document needs on first run and caches them):

```bash
# once: drop the executable into the environment you already use
python -c "import urllib.request,zipfile,io;   u='https://github.com/tectonic-typesetting/tectonic/releases/download/tectonic%400.17.0/tectonic-0.17.0-x86_64-pc-windows-msvc.zip';   zipfile.ZipFile(io.BytesIO(urllib.request.urlopen(u).read())).extractall('<your-env>/Scripts')"

# then, any time:
cd docs/report
tectonic -X compile main.tex --outdir build
```

Tectonic runs the engine as many times as the cross-references need, so one invocation is enough.
It writes nothing but `main.pdf` -- no `.aux`, no `.log` -- so building in place is clean:

```bash
tectonic -X compile main.tex --outdir .      # refreshes the checked-in main.pdf
```

**Rebuild and commit `main.pdf` whenever `main.tex`, `flowchart.tex` or a figure changes**,
otherwise the checked-in PDF quietly drifts from the source. It is ~7 MB, most of it the
figures, so each committed revision adds that much to the repository history.

**Overleaf** also works: upload this folder and compile with pdfLaTeX. Only standard TeX Live
packages are used: `graphicx, booktabs, amsmath, siunitx, xcolor, caption, subcaption, enumitem,
microtype, array, longtable, mdframed, tikz, hyperref, geometry, lmodern`.

Two preamble details worth not undoing:

- **`\DeclareSIPrefix{\micro}{\ensuremath{\mu}}`.** siunitx 3.0.x declares the micro prefix as a
  literal U+00B5, and under `T1` that slot holds *t-with-cedilla* -- so every `\micro` printed as a
  stray letter. The math mu renders correctly on every engine.
- **No `inputenc`.** It is a no-op on modern pdfLaTeX and harmful under XeTeX/LuaTeX. The source is
  pure ASCII, so nothing needs it.

## Regenerating the figures

Four of them — `calibration`, `waveform_metrics`, `preamp_clipping`, `axial_deration` — come
straight from the bundled captures:

```matlab
misetup
MakeReportFigures                              % writes to results/
MakeReportFigures([], 'docs/report/figures')   % or straight into the report
```

`MakeReportFigures` also produces a fifth figure, `reflection_check`, which is **not** included in
the document. It is the tank-reflection diagnostic the protocol points at in its Reflections
section: run it once per tank geometry to decide whether an absorber is actually needed.

The remaining figures are outputs of the wider analysis: `L11_5_verasonics_comparison` from
`PlotMI_L11_5`, and `push_axial_preamp`, `push_voltage_preamp`, `push_voltage_elements`,
`push_intensity_elements` from `PreampComparison`.

Two figures — `probe_compare_summary` and `probe_compare_bmode`, used by Appendix A — come from
outside this repository, from the `shearWaveProcessing` probe-comparison analysis
(`scripts/task3_probe_compare.py`). They are checked in as PNGs because the argument they support is
about phantom shear-wave data, not about the hydrophone. Everything else in the document is
reproducible from `data/`.

`readHWS` used to emit an HDF5 warning about a 16-byte integer type on every capture. That is
the NI 128-bit timestamp attributes, which MATLAB has no type for and which nothing here
reads; `readHWS` now silences that one identifier around its `h5info` call and restores the
caller's warning state. Do not suppress warnings globally to hide it -- other warnings from
`readHWS` (RIS sampling, a record shorter than the pulse) are real problems.

## Notes on the document

- `figures/MI_water_tank.png` was rotated 90° so it displays upright; the content is untouched.
- §2 deliberately keeps the **ideal** measurement — a full 3-D raster driven by the OWIsoft meander
  function, with the Verasonics, the stage and the digitizer synchronised — alongside the reduced
  axial-line measurement that is actually described. If that synchronisation ever gets built, §2.1
  is the specification to build against.
- Appendix A argues from *phantom shear-wave* data that the two S5-1 probes are equivalent. A
  single hydrophone capture on the MUMC probe at the known peak location would make that argument
  direct, and is worth doing.
