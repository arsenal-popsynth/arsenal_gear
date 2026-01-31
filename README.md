# arsenal gear
[![Unit Tests](https://github.com/arsenal-popsynth/arsenal_gear/actions/workflows/pytest.yml/badge.svg)](https://github.com/arsenal-popsynth/arsenal_gear/actions/workflows/pytest.yml)
[![Documentation Build](https://github.com/arsenal-popsynth/arsenal_gear/actions/workflows/documentation.yml/badge.svg)](https://github.com/arsenal-popsynth/arsenal_gear/actions/workflows/documentation.yml)
[![Linter](https://github.com/arsenal-popsynth/arsenal_gear/actions/workflows/pylint.yml/badge.svg)](https://github.com/arsenal-popsynth/arsenal_gear/actions/workflows/pylint.yml)

 Your source for Population Synthesis with an emphasis on the parameters relevant for stellar feedback from massive stars.

Documentation can be found at [arsenal-popsynth.github.io/arsenal_gear](https://arsenal-popsynth.github.io/arsenal_gear)

## Sukhbold et al. (2016) Yields & Explosion Energies

The `Sukhbold2016` yield model relies on supplementary data from:

**Sukhbold, T., Ertl, T., Woosley, S. E., Brown, J. M., & Janka, H.-T. (2016)**
"Core-Collapse Supernovae from 9 to 120 Solar Masses Based on Neutrino-powered Explosions"
ApJ 821, 38 (doi:10.3847/0004-637X/821/1/38)

**Data source**: The Astrophysical Journal Article Data
https://iopscience.iop.org/article/10.3847/0004-637X/821/1/38/meta

Download the following packages (tar.gz files):
- **Nucleosynthesis yields** (~179 MB) — contains per-progenitor `.yield_table` files with ejecta + wind contributions (stable + radioactive isotopes)
- **Explosion results** (summary files like results_N20, results_Z9.6, etc.) — includes explosion energies (E_exp in foe), fallback masses, mass cuts, etc.

After downloading and unpacking:
1. Locate the engine directories (e.g., `nucleosynthesis_yields/N20/`, `explosion_results_PHOTB/results_N20`)
2. Run the parser to generate CSVs:
