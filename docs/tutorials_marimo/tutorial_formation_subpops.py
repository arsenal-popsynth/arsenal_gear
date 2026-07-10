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
    # Case 1: Formation of stellar populations

    Assuming a fixed IMF: Salpeter between 0.08 $M_\odot$ and 100 $M_\odot$
    """
    )
    return


@app.cell(hide_code=True)
def _():
    mo.md(
        r"""
    **Choose the number of stellar populations to form**
    """
    )
    return


@app.cell
def _():
    # number of stellar populations to generate
    slider_number_subpops = mo.ui.slider(
        start=1,
        stop=3,
        step=1,
        label="Number of stellar populations to generate",
        value=1,
    )
    mo.md(f"Choose a value: {slider_number_subpops}")
    return (slider_number_subpops,)


@app.cell
def _(slider_number_subpops):
    _dict_values = {}
    for _i in range(slider_number_subpops.value):
        # create the widget to choose the initial mass of the population - solar masses
        text_mass = mo.ui.text(
            placeholder="Mass in solar masses", label="Initial mass [MSun]", value="1e5"
        )
        _dict_values[f"pop{_i+1}"] = text_mass

    dict_ui_form_props = mo.ui.dictionary(elements=_dict_values)

    mo.md(f"{dict_ui_form_props}")
    return (dict_ui_form_props,)


@app.cell(hide_code=True)
def _():
    mo.md(
        r"""
    ## Define the dictionaries with the properties of the stellar populations
    """
    )
    return


@app.cell
def _(dict_ui_form_props):
    # assume a fixed IMF - Salpeter
    _imf_min_mass = u.Quantity(0.08, unit=u.Msun)
    _imf_max_mass = u.Quantity(100, unit=u.Msun)

    dict_pops = {}
    for _key in dict_ui_form_props.value.keys():
        # SSP of single stars with a Salpeter IMF, no rotation, and no metallicity
        dict_pops[_key] = {
            "type": "single",
            "Mtot": float(dict_ui_form_props.value[_key]) * u.Msun,
            "metallicity": 0.0 * u.dimensionless_unscaled,
            "rotation": 0.0 * u.dimensionless_unscaled,
            "imf": "salpeter",
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
