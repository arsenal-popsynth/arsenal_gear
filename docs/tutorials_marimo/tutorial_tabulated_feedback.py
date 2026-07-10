# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "astropy==7.2.0",
#     "marimo",
#     "matplotlib==3.10.8",
#     "numpy==2.4.2",
#     "pandas==3.0.0",
#     "requests==2.32.5",
#     "scipy==1.17.0",
# ]
# ///

import marimo

__generated_with = "0.23.9"
app = marimo.App(
    width="full",
    app_title="ARSENAL - Formation of stellar populations",
)

with app.setup:

    # these lines work locally with uvx marimo edit, but not with export html-wasm
    import sys
    from pathlib import Path

    root = Path(__file__).resolve().parents[2]
    sys.path.insert(0, str(root))

    import astropy.units as u
    import marimo as mo  # pylint: disable=reimported
    import matplotlib.pyplot as plt
    import numpy as np

    import arsenal_gear as ag

    # load the style file for matplotlib
    plt.style.use("petroff10")
    plt.style.use(Path(__file__).parent.absolute() / "mplstyle.txt")

    import functions_plotting as fp


@app.cell(hide_code=True)
def _():
    mo.md(
        r"""
    # Case 3: Produce a table with the feedback properties for a given stellar population


    Assuming a fixed IMF: Salpeter between 0.08 $M_\odot$ and 100 $M_\odot$

    Assuming a fix stellar population (for now): 1e5MSun at 0 metallicity and null rotation
    """
    )
    return


@app.cell(hide_code=True)
def _():
    mo.md(
        r"""
    ## Define the dictionaries with the properties of the stellar populations
    """
    )
    return


@app.cell
def _():
    # assume a fixed IMF - Salpeter
    _imf_min_mass = u.Quantity(0.08, unit=u.Msun)
    _imf_max_mass = u.Quantity(100, unit=u.Msun)

    # define the time array
    _ls_time = np.logspace(6, 10, 100) * u.yr
    # define the metallicity array
    _ls_zh = [-4, -3, -2, -1, 0]

    dict_pops = {}
    for _i, _zh in enumerate(_ls_zh):
        _key = f"pop{_i}"
        # SSP of single stars with a Salpeter IMF, no rotation, and no metallicity
        dict_pops[_key] = {
            "type": "single",
            "Mtot": 1e5 * u.Msun,
            "metallicity": _zh * u.dimensionless_unscaled,
            "rotation": 0.0 * u.dimensionless_unscaled,
            "imf": "Salpeter",
            "mmin": _imf_min_mass.value,
            "mmax": _imf_max_mass.value,
            "discrete": True,
            "seed": 42,  # * int(_key.split("pop")[-1]),
            "tform": float(_key.split("pop")[-1]) * u.Myr,
        }

    # create individual populations for each SSP declared
    dict_sps = {}
    for _key in dict_pops:
        dict_sps[_key] = ag.SynthPop(pop1=dict_pops[_key], interp_op="iso")

    # collect the specific SNII rate
    dict_results = {"ndotsnii": {}}
    dict_results["time"] = _ls_time.copy()
    for _key in dict_pops:
        dict_results["ndotsnii"][_key] = (
            dict_sps[_key].ndotsn(_ls_time) / dict_sps[_key].Mtot.value
        )  # units: MSun^-1 Myr^-1
    return dict_results, dict_sps


@app.cell(hide_code=True)
def _():
    mo.md(
        r"""
    ## Create the stellar populations
    """
    )
    return


@app.cell
def _():
    return


@app.cell(hide_code=True)
def _():
    mo.md(
        r"""
    ## Plot the specifc SN rate for the population at a sample of different times
    """
    )
    return


@app.cell
def _(dict_sps, sp):
    _fig = fp.plot_vstime_sne_rate(sp, dict_sps)
    _fig.gca()
    return


@app.cell
def _(dict_sp_subpops, sp):
    _tlin = np.logspace(6, 10, 100) * u.yr

    _fig, _ax = plt.subplots(figsize=(8.4, 5))
    # plot the specific rate of SNe at the different times
    _ax.plot(_tlin, sp.ndotsn(_tlin) / sp.Mtot.value, label="Composite", c="k")
    for _key in dict_sp_subpops.keys():
        _i = int(_key.split("pop")[-1])
        _ax.plot(
            _tlin,
            dict_sp_subpops[_key].ndotsn(_tlin) / dict_sp_subpops[_key].Mtot.value,
            label=f"Pop{_i}",
            c=f"C{_i}",
        )

    _fig.gca()
    return


@app.cell(hide_code=True)
def _():
    mo.md(
        r"""
    **(3) Compute a history of the population's bolometric luminosity for a list of times.**
    """
    )
    return


@app.cell
def _(dict_sp_subpops, sp):
    _fig = fp.plot_vstime_bolometric_luminosity(sp, dict_sp_subpops)
    _fig.gca()
    return


@app.cell(hide_code=True)
def _():
    mo.md(
        r"""
    **(4) Compute a history of the population's effective temperature for a list of times.**
    """
    )
    return


@app.cell
def _(dict_sp_subpops, sp):
    _fig = fp.plot_vstime_effective_temperature(sp, dict_sp_subpops)
    _fig.gca()
    return


@app.cell(hide_code=True)
def _():
    mo.md(
        r"""
    # MISSING: yields! winds!
    """
    )
    return


if __name__ == "__main__":
    app.run()
