[README.md](https://github.com/user-attachments/files/28859362/README.md)
# QSO Continuum Fitting Code

Code used in **Mishra & Muzahid (2022)**, *The Astrophysical Journal*.

This repository provides the continuum fitting pipeline for quasar (QSO) spectra, designed to work with SDSS spectral data products. It fits a spline-based continuum to the MgII spectral region and identifies whether known emission lines fall within the fitting window.

---

## Reference

If you use this code, please cite:

> Mishra, S. & Muzahid, S. (2022). *Title of Paper*. The Astrophysical Journal.  
> DOI: [add DOI here]

---

## Overview

The pipeline:

1. Opens and reads calibrated FITS spectra (wavelength, flux, error, and pre-existing continuum estimates).
2. Defines a fitting window around the MgII λλ2796, 2803 doublet in the observed frame, covering a velocity range of ±30,000 km/s.
3. Computes the signal-to-noise ratio (SNR) in the MgII window.
4. Checks whether major emission lines (Lyα, Si IV, C IV, C III], MgII) from the quasar fall within the fitting window.
5. Applies sigma clipping and a weighted cubic spline fit to model the continuum.
6. Writes per-spectrum diagnostics to an output log file.

---

## Repository Structure

```
.
├── continuum_fit.py       # Main fitting script
├── main_1-5r500_sample-4-NoGap-NoBAL.txt   # Input pair catalog (not included; provide your own)
└── README.md
```

---

## Dependencies

Install all required packages via `pip`:

```bash
pip install numpy scipy astropy scikit-learn matplotlib lmfit
```

| Package        | Purpose                                      |
|----------------|----------------------------------------------|
| `numpy`        | Array operations                             |
| `scipy`        | Spline fitting, signal processing, optimization |
| `astropy`      | FITS I/O, sigma clipping                     |
| `scikit-learn` | Linear model utilities                       |
| `matplotlib`   | Plotting and PDF output                      |
| `lmfit`        | Parameter fitting interface                  |

---

## Input

### FITS Spectra

Each spectrum must be a FITS file with a binary table extension (HDU index 1) containing the following columns:

| Column         | Description                          |
|----------------|--------------------------------------|
| `WAVE`         | Observed wavelength array (Å)        |
| `FLUX`         | Flux array                           |
| `ERROR`        | 1σ flux error array                  |
| `CONTI_SDSS`   | Pre-existing SDSS continuum estimate |
| `CONTI_SPLINE` | Pre-existing spline continuum        |

The spectra directory path is set in the script via the `path` variable:

```python
path = '/path/to/your/spectra/'
```

### Pair Catalog

A whitespace-delimited text file (with a one-line header) listing QSO–cluster pairs. The relevant columns read by the script are:

| Column index | Variable      | Description                        |
|-------------|---------------|------------------------------------|
| 0           | `ra_q`        | QSO right ascension                |
| 1           | `dec_q`       | QSO declination                    |
| 2           | `zemi_q`      | QSO emission redshift              |
| 4           | `SNR`         | Signal-to-noise ratio              |
| 5           | `Mi`          | Absolute i-band magnitude          |
| 6           | `imag`        | Apparent i-band magnitude          |
| 7           | `file_in`     | FITS filename                      |
| 8           | `ra_cl`       | Cluster right ascension            |
| 9           | `dec_cl`      | Cluster declination                |
| 10          | `z_cl`        | Cluster redshift                   |
| 11          | `M500`        | Cluster mass M500                  |
| 12          | `R500`        | Cluster radius R500                |
| 13          | `q_cl_vel`    | QSO–cluster velocity separation    |
| 14          | `rhy`         | Projected separation (Mpc)         |
| 15          | `rhy_R500`    | Projected separation in R500 units |
| 17          | `RL500`       | Richness within projected R500     |
| 18          | `N500sp`      | Spectroscopic N500                 |
| 19          | `N500`        | Total N500                         |
| 20          | `R500_wen`    | R500 from Wen catalog              |
| 22          | `pairid`      | Pair ID (int)                      |
| 23          | `gal_flg_spec`| Galaxy spectroscopic flag (int)    |
| 24          | `gal_flghot`  | Galaxy hot-gas flag (int)          |
| 3           | `name_q`      | QSO name                           |
| 16          | `cat`         | Catalog name                       |
| 21          | `name_cl`     | Cluster name                       |

The default catalog filename is set as:

```python
pair_file = 'main_1-5r500_sample-4-NoGap-NoBAL.txt'
```

---

## Output

The script appends results to a plain-text log file called `tmp4_em` with the following columns:

```
File   pix   snr_er   snr_f   lam_low   lam_up   em-info
```

| Column    | Description                                                         |
|-----------|---------------------------------------------------------------------|
| `File`    | FITS filename                                                        |
| `pix`     | Number of pixels in the MgII fitting window                         |
| `snr_er`  | Median SNR using flux error array                                   |
| `snr_f`   | Median SNR estimated from flux standard deviation                   |
| `lam_low` | Lower wavelength boundary of fitting window (Å)                     |
| `lam_up`  | Upper wavelength boundary of fitting window (Å)                     |
| `em-info` | Whether a major QSO emission line falls within the window (`Emission-In` / `Emission-Out`) |

---

## Key Functions

### `open_calibrate_fits(filename, path)`
Opens a calibrated FITS file and returns the wavelength, flux, error, and continuum arrays.

### `spline_fit(ts, ys, ys_err, knot)`
Fits a weighted cubic B-spline to the spectrum using interior knots placed at quantile positions. Weights are set to `1/σ²`.

### `sigma_clip_box(datx, nsig_up, nsig_low, box)`
Performs iterative sigma clipping within equal-width boxes along the spectrum. Returns indices of unclipped data points.

### `running_median(datx, daty, daty_err, bin_size=21)`
Computes a running weighted average and unweighted median over a sliding window.

### `main_fit_part(runid, file_in, zemi_q, z_cl, path)`
Main routine that processes a single QSO spectrum: loads data, defines the fitting window, checks for emission line contamination, computes SNR, and writes results to the output log.

---

## Usage

Edit the `path` and `pair_file` variables at the bottom of the script to point to your data, then run:

```bash
python continuum_fit.py
```

Results are appended to `tmp4_em` in the current working directory.

---

## Notes

- The MgII fitting window spans ±30,000 km/s around MgII λ2796.35 at the cluster redshift, anchored between the quasar redshift lower bound and the cluster upper bound.
- Emission line contamination is checked for Lyα (1215.67 Å), Si IV (1398.27 Å), C IV (1549.49 Å), C III] (1908.73 Å), and MgII (2799.94 Å), with a proximity threshold of ±3000 km/s.
- The code is designed for SDSS DR spectra processed through a PCA+spline continuum pipeline as described in Mishra & Muzahid (2022).

---

## Contact

**Sapna Mishra**  
[Add institutional email or GitHub profile link here]
