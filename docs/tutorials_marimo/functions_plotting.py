### Routines to show properties of the stellar populations
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np

import arsenal_gear as ag


def plot_vstime_effective_temperature(sp, dict_sp_subpops):
    """Routine to show the time evolution of the effective temperature of a composite stellar population, and its individual populations
    Input:
    sp: composite stellar population
    dict_sp_subpops: dictionary containing the individual stellar populations
    Output:
    _fig: figure
    """
    # define the range of times
    _tlin = np.logspace(5.1, 9, 20) * u.yr

    _fig, _ax = plt.subplots(figsize=(8.4, 5))
    _ax.plot(_tlin, sp.teff(_tlin), "o-", label="Composite", c="k")
    for _key in dict_sp_subpops.keys():
        _i = int(_key.split("pop")[-1])
        _ax.plot(_tlin, dict_sp_subpops[_key].teff(_tlin), label=f"Pop{_i}", c=f"C{_i}")

    # format the axes
    _ax.legend()
    _ax.set_yscale("log")
    _ax.set_xscale("log")
    _ax.set_ylabel(r"$T_{\rm eff} \, [K]$")
    _ax.set_xlabel("Time [yr]")
    return _fig


def plot_vstime_bolometric_luminosity(sp, dict_sp_subpops):
    """Routine to show the time evolution of the bolometric luminosity of a composite stellar population, and its individual populations
    Input:
    sp: composite stellar population
    dict_sp_subpops: dictionary containing the individual stellar populations
    Output:
    _fig: figure
    """
    # define the range of times
    _tlin = np.logspace(5.1, 9, 20) * u.yr

    _fig, _ax = plt.subplots(figsize=(8.4, 5))
    _ax.plot(_tlin, sp.lbol(_tlin), "o-", label="Composite", c="k")
    for _key in dict_sp_subpops.keys():
        _i = int(_key.split("pop")[-1])
        _ax.plot(_tlin, dict_sp_subpops[_key].lbol(_tlin), label=f"Pop{_i}", c=f"C{_i}")

    # format the axes
    _ax.legend()
    _ax.set_yscale("log")
    _ax.set_xscale("log")
    _ax.set_ylabel(r"$L_{\rm bol} \, [L_{\odot}]$")
    _ax.set_xlabel("Time [yr]")
    return _fig


def plot_vstime_sne_rate(sp, dict_sp_subpops):
    """Routine to show the time evolution of the SNe rate of a composite stellar population, and its individual populations
    Input:
    sp: composite stellar population
    dict_sp_subpops: dictionary containing the individual stellar populations
    Output:
    _fig: figure
    """
    # define an array of times to compute the specific rate of SNe at
    _tlin = np.logspace(6, 8, 100) * u.yr

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

    # format the axes
    _ax.legend()
    _ax.set_xscale("log")
    _ax.set_yscale("log")
    _ax.set_xlabel("Time [yr]")
    _ax.set_ylabel(r"$\dot{N}_{\rm SN}\, [M_{\odot}^{-1}\,\,\mathrm{Myr}^{-1}]$")
    _ax.set_xlim(1e6, 1e8)
    _ax.set_ylim(1e-5, 8e-4)
    return _fig


def plot_combined_stochastic_imfs(dict_sp_subpops):
    """Routine to show the stochastic IMFs sampled within a composite stellar population, and its stellar populations (as long as they are discrete)
    Input:
    sp: composite stellar population
    dict_sp_subpops: dictionary containing the individual stellar populations
    Output:
    _fig: figure
    """

    def add_stochastic_imf(
        fig: plt.Figure,
        ax: plt.Axes,
        ind: int,
        masses: np.ndarray,
        min_mass=0.08 * u.Msun,
        max_mass=100 * u.Msun,
        imf=ag.formation.dist_funcs.imf,
    ):
        """Routine to add the sampled IMF of a stellar population in the combined figure
        Input:
        fig: Figure object
        ax: axes object
        masses: discrete masses
        min_mass: minimum mass sampled
        max_mass: maximum mass sampled
        imf: IMF from which the mass has been sampled
        Output:
        _fig: figure"""
        # make a histogram of the sampled masses and compare to the underlying PDF

        # create the bins for the histogram
        mbins = np.logspace(np.log10(min_mass.value), np.log10(max_mass.value), 50)
        # create instance for the figure
        # add the histogram of the sampled masses and the underlying PDF
        ax.hist(
            masses,
            bins=mbins,
            density=True,
            histtype="step",
            color=f"C{ind}",
            label=f"Sampled Pop{ind}",
        )
        ax.plot(mbins, imf.pdf(mbins), color="C3")

        return fig

    fig, ax = plt.subplots(figsize=(8.4, 5))
    ls_pops = list(dict_sp_subpops.keys())
    ls_pops.sort()
    for i, key in enumerate(ls_pops):
        # create the same IMF as used to create the SSP
        _imf = ag.formation.dist_funcs.imf.IMF.get_imf(
            imf_identifier=dict_sp_subpops[key].fdict["pop1"].get("imf", "salpeter"),
            min_mass=dict_sp_subpops[key].fdict["pop1"].get("mmin", 0.08) * u.Msun,
            max_mass=dict_sp_subpops[key].fdict["pop1"].get("mmax", 100) * u.Msun,
            seed=dict_sp_subpops[key].fdict["pop1"].get("seed", 42),
        )
        _masses = dict_sp_subpops[key].discrete_masses.value
        add_stochastic_imf(
            fig,
            ax,
            ind=i,
            masses=_masses,
            min_mass=dict_sp_subpops[key].fdict["pop1"].get("mmin", 0.08) * u.Msun,
            max_mass=dict_sp_subpops[key].fdict["pop1"].get("mmax", 100) * u.Msun,
            imf=_imf,
        )

        # format the axes
        ax.set_yscale("log")
        ax.set_xscale("log")
        ax.set_xlim(_imf.min_mass, _imf.max_mass)
        ax.set_ylim(1e-6, 1e2)
        ax.set_xlabel(r"$M \, [M_{\odot}]$")
        ax.set_ylabel(r"$dN/d\log M$")
        ax.legend(loc="best")

    return fig


def plot_stochastic_imf(
    masses: np.ndarray,
    min_mass=0.08 * u.Msun,
    max_mass=100 * u.Msun,
    imf=ag.formation.dist_funcs.imf,
):
    """Routine to show the stochastic IMFs sampled within a discrete stellar population
    Input:
    masses: discrete masses
    min_mass: minimum mass sampled
    max_mass: maximum mass sampled
    imf: IMF from which the mass has been sampled
    Output:
    _fig: figure
    """
    # make a histogram of the sampled masses and compare to the underlying PDF

    # create the bins for the histogram
    mbins = np.logspace(np.log10(min_mass.value), np.log10(max_mass.value), 50)
    # create instance for the figure
    fig, ax = plt.subplots(figsize=(8.4, 5))
    # add the histogram of the sampled masses and the underlying PDF
    ax.hist(
        masses,
        bins=mbins,
        density=True,
        histtype="step",
        color="C0",
        label="Sampled IMF",
    )
    ax.plot(mbins, imf.pdf(mbins), color="C3", label="IMF PDF")
    # format the axes
    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.set_xlim(min_mass.value, max_mass.value)
    ax.set_ylim(1e-6, 1e2)
    ax.set_xlabel(r"$M \, [M_{\odot}]$")
    ax.set_ylabel(r"$dN/d\log M$")
    ax.legend(loc="best")
    return fig
