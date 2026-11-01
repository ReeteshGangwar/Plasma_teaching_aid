# Plasma Spectroscopy Teaching Tools

Two single-file browser applications for teaching grating spectrometers and optical emission spectroscopy (OES) of argon plasmas.

| Tool | Folder | What it does |
|---|---|---|
| **Czerny–Turner Bench** | [`czerny-turner-bench/`](czerny-turner-bench/) | Interactive Czerny–Turner spectrometer: 3D model, cut-away animation and to-scale ray geometry. Includes live formulas (grating equation, dispersion, resolution budget), simulated Ar–H plasma and Hg–Ar lamp spectra, and Doppler/Stark line-width exercises. |
| **Plasma OES Analyzer** | [`oes-analyzer/`](oes-analyzer/) | Analyses student OES data: triplicate spectra with background subtraction, peak identification, an Ar I Boltzmann plot with LTE checks, N₂ rotational temperature, and electron density from Ar I 696.5 nm broadening. |

Everything runs in the browser. No server, no installation, and data files are never uploaded.

## Using the tools

**Online (recommended):** enable GitHub Pages for this repository (*Settings → Pages → Deploy from branch → `main` / root*). Students then open:

- `https://<your-user>.github.io/<repo>/` (landing page)
- `https://<your-user>.github.io/<repo>/czerny-turner-bench/`
- `https://<your-user>.github.io/<repo>/oes-analyzer/`

**Offline:** download `oes-analyzer/index.html` and double-click it. It contains its plotting library and needs no internet. The Czerny–Turner Bench needs internet the first time it loads, to fetch three.js and MathJax from public CDNs.

## Plasma OES Analyzer

### Input data

- **File format:** plain text, two columns: wavelength (nm) and intensity. Header lines are skipped, and comma, tab, semicolon or space separators all work, as do decimal commas. Only the first two numeric columns are used.
- **Per condition:**
  - 3 plasma spectra and 3 background spectra (dark plus ambient light, plasma off).
  - Optionally, a high-resolution record of the Ar I 696.5 nm line (plasma and background) for the electron-density tab.
- **Automatic sorting:** if all files of one condition are dropped together, they are sorted by filename. Names containing `bg`, `background`, `dark` or `blank` go to background; `696`, `hires` or `hr` go to the high-resolution set; everything else goes to plasma.
- **Intensity units:** spectra are assumed to be already corrected for spectral response (energy units). The *photons* setting exists for uncorrected counts.
- **Example files:** [`oes-analyzer/examples/`](oes-analyzer/examples/) has synthetic example files in the expected format. The same data load with the *Load example data* button.

### Analysis steps (tabs)

1. **Data:**
   - Net spectrum = each plasma replicate − mean background, scaled by exposure time if both exposures are given.
   - The result is reported as mean ± SD over the replicates.
   - Saturated pixels are flagged.
2. **Spectrum & peaks:**
   - Automatic wavelength-shift correction from strong Ar I lines.
   - Peak finding with the threshold set from the background noise.
   - Peaks are assigned to Ar I/II, Hα/Hβ/Hγ, O I, OH A–X, N₂ second positive system (SPS) and N₂⁺ first negative system (FNS).
3. **Boltzmann & LTE check:**
   - ln(Iλ/gA) against Eₖ for Ar I lines up to 900 nm (5p→4s, 4p→4s, 5d/6d/7s→4p).
   - Three intensity sources: integrated area, peak-finder area or height, or manual entry (typed or pasted).
   - LTE checks: straight-line fit quality, 4p-only against higher-level slopes, lines sharing the same upper level, and the McWhirter criterion using nₑ from tab 5.
   - When the checks fail, the tab explains why a collisional–radiative (CR) model is needed.
4. **Gas temperature (N₂):** fit of the N₂ C³Πᵤ→B³Π_g (0,0) band at 337 nm, giving T_rot per replicate, a Δχ² uncertainty and residuals. There is also an export for [massiveOES](https://ics.muni.cz/en/research/publications/1437096).
5. **nₑ from Ar 696.5 nm:**
   - Voigt fit with the Gaussian part fixed from the instrument width and Doppler width.
   - W_L = W_Stark + W_vdW, which gives nₑ.
   - nₑ is reported as an upper limit when the Stark width is below its uncertainty.
6. **Results:** all conditions side by side, plus intensity ratios relative to Ar I 750.4 nm (OH, N₂, N₂⁺, Hβ, Hα, O I). Exports to CSV.

### Scientific assumptions: please read

- **Ar I atomic data:**
  - Lines tagged *NIST* were taken from the NIST Atomic Spectra Database.
  - Lines tagged *verify A,g,E* should be checked against [NIST ASD](https://physics.nist.gov/PhysRefData/ASD/lines_form.html) before quantitative use.
  - The 545–565 nm lines follow Ohno et al., *Plasma Fusion Res.* **1**, 028 (2006).
  - All constants are editable in the app.
- **N₂ rotational model:**
  - Simplified Π–Π model with unresolved triplet structure; B′ = 1.8154 cm⁻¹, B″ = 1.6285 cm⁻¹, Hönl–London factors, and a Gaussian instrument function.
  - Suitable for teaching and for comparing conditions measured on the same instrument.
  - Use massiveOES for publication-quality temperatures.
  - T_rot(N₂) can overestimate the gas temperature in argon because of energy transfer from Ar metastables.
- **Ar I 696.5 nm Stark coefficient:**
  - The default, w_S = 0.0112 nm per 10¹⁶ cm⁻³ at Tₑ ≈ 10⁴ K, reproduces the conversion of Lu et al., *Sci. Rep.* **4**, 7488 (2014).
  - Check it against Konjević & Konjević, *Spectrochim. Acta B* **52**, 2077 (1997).
- **Ar I 696.5 nm van der Waals width:**
  - C = 7.5×10⁻⁷ λ² ⟨R²⟩^(2/5) (Yubero et al.; see *Atoms* **5**, 41, 2017), with ⟨R²⟩ ≈ 36.6 a₀² from the Coulomb approximation.
  - At 1 atm this width (~0.025 nm) dominates the Stark width below ~10¹⁵ cm⁻³.
- **Excitation temperature:** the Boltzmann plot gives T_exc, which equals Tₑ only under (partial) LTE. Atmospheric-pressure argon jets usually fail the McWhirter criterion, in which case use a CR model, e.g.:
  - X.-M. Zhu & Y.-K. Pu, *J. Phys. D* **43**, 015204 (2010)
  - S. Siepa et al., *J. Phys. D* **47**, 445201 (2014)

### Rebuilding the analyzer from source

`oes-analyzer/index.html` is generated from the files in `oes-analyzer/src/`:

| File | Contents |
|---|---|
| `engine.js` | Parsing, processing, peak finding, Boltzmann plot, LTE checks, N₂ model and fits, Voigt and broadening, synthetic example data |
| `app.js` | User interface |
| `ui.html` | Page structure and styles |
| `build.py` | Inlines Plotly and writes `../index.html` |

```bash
cd oes-analyzer/src
npm pack plotly.js-dist-min@2.27.0 && tar xzf plotly.js-dist-min-2.27.0.tgz
python3 build.py
```

## Czerny–Turner Bench

Open `czerny-turner-bench/index.html`. It has three views:

- **3D model:** drag to orbit, scroll or pinch to zoom. Includes an exploded view.
- **Cut-away:** the animated light path in 2D.
- **To-scale:** the true ray geometry for the current settings.

Settings on the left (grating, order, focal length, f-number, deviation angle, slit, detector, light source) update the formulas, the resolution budget, the simulated CCD spectrum and the line close-up.

## Licence

Code: MIT (see [`LICENSE`](LICENSE)). Third-party libraries keep their own licences (see [`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md)).
