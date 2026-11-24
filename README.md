# Variable Star Searcher

*Find uncatalogued variable stars in TESS, Kepler and K2 light curves, on a home server*

A companion of Exoplanet Searcher that runs as a Docker service at the lowest priority: it searches TESS (2-min and full-frame), Kepler and K2 light curves for variable stars that are **not yet catalogued**, checks them against VSX and other catalogues, and prepares the data needed for a submission to the AAVSO VSX.

## Features

* One thread and at most 3 GB of RAM; FITS files deleted after each star
* Periodogram and classification, blending check, epoch in HJD
* VSX-ready names (HIP/TYC/2MASS/TIC) and Gaia coordinates via TIC
* 1.2: HEASARC queries try CSV, then VOTable (the old `votable/td` format was dropped by HEASARC)

> **Language note.** The interface is in **Italian** (large, adjustable text); plots, PDF/PNG cards and scientific output are in **English**. Code comments are in Italian. The original Italian guide is [`LEGGIMI.txt`](LEGGIMI.txt).

**Version:** 1.2 (2026-11-07). **Status:** 1.2 written; 1.1 ran on a NAS until HEASARC changed formats.

## Quick start

Requirements: Python 3.10+ (tested with Python 3.14 on Ubuntu) and Docker for the server part.

```bash
git clone https://github.com/EDOKILLER/variable-star-searcher.git
cd variable-star-searcher
python3 -m pip install -r requirements.txt

cd nas
python3 servizio_var.py --prova     # offline self-test
# then deploy nas/docker-compose.yml (see nas/LEGGIMI_NAS.txt)
```

### Server / NAS part

The `nas/` folder runs in a `python:3.12` Docker container (tested on a UGREEN NAS with UGOS Docker, any Linux box works). Copy the files into the shared folder, edit the volume path in `nas/docker-compose.yml`, and deploy it; the container installs its libraries on first start and runs at the lowest CPU and disk priority. Step-by-step Italian instructions: `nas/LEGGIMI_NAS.txt`.

## Data sources and acknowledgements

Data are downloaded at run time from public services. If you publish results, please follow their acknowledgement policies:

* [Mikulski Archive for Space Telescopes (STScI): Hubble, JWST, Kepler, TESS, GALEX, Pan-STARRS](https://archive.stsci.edu/)
* [NASA High Energy Astrophysics Science Archive Research Center](https://heasarc.gsfc.nasa.gov/)
* [VizieR catalogue service, CDS, Strasbourg (the Harvard mirror is used)](https://vizier.cds.unistra.fr/)
* [ESA Gaia archive](https://gea.esac.esa.int/archive/)

This project is independent and not endorsed by any of these institutions.

## License

Code: [BSD 3-Clause](LICENSE) © 2026 Edoardo Battistel. You can use, modify and redistribute it, also commercially, as long as you keep the copyright notice; the author's name may not be used to promote derived works without permission. If you use it in research, see [`CITATION.cff`](CITATION.cff).

---
Part of a series of astronomy programs by Edoardo Battistel — see also [X-rays from Black Holes](https://github.com/EDOKILLER/x-rays-from-black-holes) and the other repositories on [this profile](https://github.com/EDOKILLER?tab=repositories).
