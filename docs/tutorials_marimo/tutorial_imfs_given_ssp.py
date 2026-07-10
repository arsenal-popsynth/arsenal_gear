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

    import arsenal_gear as ag

    # load the style file for matplotlib
    plt.style.use("petroff10")
    plt.style.use(Path(__file__).parent.absolute() / "mplstyle.txt")

    import functions_plotting as fp


@app.cell(hide_code=True)
def _():
    mo.md(
        r"""
    # Case 2: Creating stellar populations with different IMFs

    Assuming a fixed stellar population: a $10^5~M_\odot$ population of single stars, with no metallicity nor rotation
    """
    )
    return


@app.cell
def _():
    # parameters of the IMF - lowest/highest masses and slope
    text_imf_min = mo.ui.text(
        placeholder="Lowest mass IMF",
        label="Lowest mass of the IMF [MSun]",
        value="0.08",
    )
    text_imf_max = mo.ui.text(
        placeholder="Highest mass IMF",
        label="Highest mass of the IMF [MSun]",
        value="100",
    )

    return text_imf_max, text_imf_min


@app.cell
def _(text_imf_min):
    mo.md(f"{text_imf_min}")
    return


@app.cell
def _(text_imf_max):
    mo.md(f"{text_imf_max}")
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
def _(text_imf_max, text_imf_min):
    ls_imfs = ["Salpeter", "Kroupa1993", "Kroupa2001", "MillerScalo", "Chabrier"]
    # assume a fixed IMF - Salpeter
    _imf_min_mass = u.Quantity(text_imf_min.value, unit=u.Msun)
    _imf_max_mass = u.Quantity(text_imf_max.value, unit=u.Msun)

    dict_pops = {}
    for _i, _imf_name in enumerate(ls_imfs):
        _key = f"pop{_i+1}"
        # SSP of single stars with a Salpeter IMF, no rotation, and no metallicity
        dict_pops[_key] = {
            "type": "single",
            "Mtot": 1e5 * u.Msun,
            "metallicity": 0.0 * u.dimensionless_unscaled,
            "rotation": 0.0 * u.dimensionless_unscaled,
            "imf": _imf_name,
            "mmin": _imf_min_mass.value,
            "mmax": _imf_max_mass.value,
            "discrete": True,
            "seed": 42 ** int(_key.split("pop")[-1]),
            "tform": float(_key.split("pop")[-1]) * u.Myr,
        }
    return (dict_pops,)


@app.cell(hide_code=True)
def _():
    mo.md(
        r"""
    **Create the stellar populations**
    """
    )
    return


@app.cell
def _(dict_pops):
    # create a composite stellar population
    sp = ag.SynthPop(**dict_pops, interp_op="iso")
    # for comparison, create individual populations for each SSP declared
    dict_sp_subpops = {}
    for _key in dict_pops:
        dict_sp_subpops[_key] = ag.SynthPop(pop1=dict_pops[_key], interp_op="iso")
    return dict_sp_subpops, sp


@app.cell(hide_code=True)
def _():
    mo.md(
        r"""
    **(1) Check the distribution of masses in comparison to the underlying PDF from which it is sampled**
    """
    )
    return


@app.cell
def _(dict_sp_subpops):
    _fig = fp.plot_combined_stochastic_imfs(dict_sp_subpops)
    _fig.gca()
    return


@app.cell(hide_code=True)
def _():
    mo.md(
        r"""
    **(2) Plot the specifc SN rate for the population at a sample of different times**
    """
    )
    return


@app.cell
def _(dict_sp_subpops, sp):
    _fig = fp.plot_vstime_sne_rate(sp, dict_sp_subpops)
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
