"""
binaries
==========

This submodule contains all the code required to sample from mass-dependent a
binary fraction and distributions of orbital parameters.
"""

import astropy.units as u
import numpy as np
from astropy.units import Quantity
from scipy.stats import loguniform, rv_continuous, uniform


class BinaryDistribution:
    """
    This class is the superclass of all distributions of binary properties
    in which the binary fraction and distributions of orbtial periods and
    mass ratios are not independent. This samples from a 3D probability
    distribution that incldues both binaries and single stars.
    This can also include a user-specified mass limit above which all
    binaries have orbital periods log(P/days) <= 3.8, to account for the
    inner companion in hierarchical systems.

    :param model: BPASS model, see BPASS manual for details
    :type model: str
    :param m_inner: Mass above which log(P/days) <= 3.8
    :type m_inner: astropy mass unit
    :param model: BPASS directory (likely on scratch, given the size)
    :type model: str

    """

    def __init__(self, model: str, m_inner: Quantity["mass"], bpass_dir: str):
        self.models = [
            "100_100",
            "100_300",
            "135_100",
            "135_300",
            "135all_100",
            "170_100",
            "170_300",
            "chab100",
            "chab300",
        ]
        if model not in self.models:
            raise ValueError(f"Model {model} does not exist.")
        if m_inner < (9 * u.Msun):
            raise Warning("Please select a value >= 9 Msun.")
        self.model: str = model
        self.m_inner: float = m_inner.to(u.Msun).value
        self.dir: str = bpass_dir + "/"

    def generate_probabilities(self):
        """
        Generate the appropriate 3D grid of probabilities from the
        BPASS input file.

        """

        def get_probability_matrix(file_path, m1_vals, q_vals, logp_vals):

            def get_prob_single(ln):
                m1 = float(ln[ln.rfind("-") + 1 :])
                save_w = np.array([np.argmin(np.abs(m1_vals - m1)), 0, -1])
                return save_w

            def get_prob_binary(ln):
                ind_p = ln.rfind("-")
                ind_q = ln[:ind_p].rfind("-")
                ind_m = ln[:ind_q].rfind("-")
                logp = float(ln[ind_p + 1 :])
                q = float(ln[ind_q + 1 : ind_p])
                m1 = float(ln[ind_m + 1 : ind_q])
                save_w = np.array(
                    [
                        np.argmin(np.abs(m1_vals - m1)),
                        np.argmin(np.abs(q_vals - q)),
                        np.argmin(np.abs(logp_vals - logp)),
                    ]
                )
                return save_w

            with open(file_path, encoding="utf-8") as f:
                prob = np.zeros((len(m1_vals), len(q_vals), len(logp_vals)))
                save_n = False
                save_w = np.zeros(3)
                for line in f:
                    ln = line.strip()
                    if "NEWSINMODS" in ln:
                        # e.g. NEWSINMODS/z014/sneplot-z014-300
                        # Read only mass
                        save_n = True
                        save_w = get_prob_single(ln)

                    elif "NEWBINMODS/NEWBINMODS" in ln:
                        # e.g. NEWBINMODS/NEWBINMODS/z014/sneplot-z014-300-0.1-0
                        # Read period, mass ratio, and mass
                        save_n = True
                        save_w = get_prob_binary(ln)

                    elif save_n:
                        save_n = False
                        ln = np.array(ln.split())
                        if prob[save_w[0], save_w[1], save_w[2]] > 0:
                            prob[save_w[0], save_w[1], save_w[2]] += float(ln[0])
                        else:
                            prob[save_w[0], save_w[1], save_w[2]] = ln[0]

            return prob

        def edit_probability_matrix(prob, m_inner, m1_vals, q_vals, logp_vals):

            # If mass limit for inner companion
            # Zero out probability above 3.8 if M > m_inner
            _mask_m = np.where(m1_vals >= m_inner)[0]
            _mask_q = np.where(q_vals > 0)[0]
            _mask_p = np.where(logp_vals > 3.8)[0]
            _mask_i = np.where(logp_vals <= 3.8)[0]
            _norm_inner = (prob[_mask_m[0] :, _mask_q[0] :, : _mask_i[-1]]).sum()
            _norm_outer = (prob[_mask_m[0] :, _mask_q[0] :, _mask_p[0] :]).sum()
            prob[_mask_m[0] :, _mask_q[0] :, _mask_p[0] :] = np.zeros(
                (len(_mask_m), len(_mask_q), len(_mask_p))
            )
            prob[_mask_m[0] :, _mask_q[0] :, : _mask_i[-1]] *= (
                _norm_inner + _norm_outer
            ) / _norm_inner

            return prob

        m1_vals = np.concatenate(
            (
                np.arange(0.1, 0.7, 0.1),
                np.arange(0.8, 2.2, 0.1),
                np.array([2.3, 2.5, 2.7, 3, 3.2, 3.5, 3.7]),
                np.arange(4, 10, 0.5),
                np.arange(10, 26, 1),
                np.array([30, 35, 40, 50, 60, 70, 80, 100, 120, 150, 200, 300]),
            )
        )

        q_vals = np.arange(0, 1, 0.1)  # q = 0 for single
        logp_vals = np.arange(0, 4.4, 0.2)  # 4.2 w/ q = 0 for single

        file_path = (
            self.dir
            + "bpass_v2.2.1_imf"
            + self.model
            + "/input_bpass_z014_bin_imf"
            + self.model
        )

        prob = get_probability_matrix(file_path, m1_vals, q_vals, logp_vals)

        if self.m_inner < int(self.model[-3:]):

            prob = edit_probability_matrix(
                prob, self.m_inner, m1_vals, q_vals, logp_vals
            )
            np.save(self.dir + "prob_bin_imf" + self.model + ".npy", prob)

        else:
            np.save(self.dir + "prob_bin_imf" + self.model + ".npy", prob)

    def discrete_sampling(
        self, mtot: Quantity["mass"]
    ) -> tuple[Quantity["mass"], Quantity["mass"], Quantity["mass"], Quantity["time"]]:
        """
        Draw individual binaries and single stars from a BPASS input file
        with a given IMF slope and upper mass cut-off, with a target total mass

        :param mtot: Target total mass of the sample
        :type mtot:  Quantity["mass"]
        :return:     Single stars' masses
        :rtype:      Quantity["mass"]
        :return:     Primary masses
        :rtype:      Quantity["mass"]
        :return:     Orbital periods
        :rtype:      Quantity["time"]

        """

        def sample_distribution(prob, mtot, m, q, logp):

            # Sample by flattening the array
            # Assume that the average system mass is below 1.5 MSun
            num_samples = int(1.5 * mtot.to(u.Msun).value)
            prob_flat = prob.ravel()
            samp_inds = np.random.choice(
                a=np.prod(prob.shape), size=num_samples, p=prob_flat
            )
            # Convert flattened indices back to 3D coordinates
            indices = np.unravel_index(samp_inds, prob.shape)
            samp_m = m[indices[0]]
            samp_q = q[indices[1]]
            samp_p = logp[indices[2]]

            return samp_m, samp_q, samp_p

        def truncate_at_target_mass(mtot, samp_m, samp_q, samp_p):

            _target_mass = np.argmin(
                np.abs(mtot.to(u.Msun).value - (samp_m * (1 + samp_q)).cumsum())
            )
            _primary_masses = samp_m[:_target_mass] * u.Msun
            _companion_masses = samp_q[:_target_mass] * samp_m[:_target_mass] * u.Msun
            _orbital_periods = 10 ** samp_p[:_target_mass] * u.d

            _mask_s = np.where(_companion_masses == 0)[0]
            _mask_b = np.where(_companion_masses > 0)[0]

            s_masses = _primary_masses[_mask_s]
            p_masses = _primary_masses[_mask_b]
            c_masses = _companion_masses[_mask_b]
            periods = _orbital_periods[_mask_b]

            return s_masses, p_masses, c_masses, periods

        prob = np.load(self.dir + "prob_bin_imf" + self.model + ".npy")
        prob = prob / prob.sum()

        m1 = np.concatenate(
            (
                np.arange(0.1, 0.7, 0.1),
                np.arange(0.8, 2.2, 0.1),
                np.array([2.3, 2.5, 2.7, 3, 3.2, 3.5, 3.7]),
                np.arange(4, 10, 0.5),
                np.arange(10, 26, 1),
                np.array([30, 35, 40, 50, 60, 70, 80, 100, 120, 150, 200, 300]),
            )
        )

        q = np.arange(0, 1, 0.1)  # q = 0 for single
        logp = np.arange(0, 4.4, 0.2)  # 4.2 w/ q = 0 for single

        samp_m, samp_q, samp_p = sample_distribution(prob, mtot, m1, q, logp)
        s_masses, p_masses, c_masses, periods = truncate_at_target_mass(
            mtot, samp_m, samp_q, samp_p
        )

        return s_masses, p_masses, c_masses, periods


class Fraction:
    """
    This class is the superclass of all mass-dependent binary fraction
    This assumes that the binary fraction is a step function

    :param fraction: Binary fraction of the mass bins, of length k
    :type fraction: float
    :param mass_bins: Limit of the mass bins, of length k-1
    :type mass_bins: astropy mass unit
    :param stars: Potential primaries
    :type stars: Quantity["mass"], list of stellar masses
    :param name: Name of the binary fraction function
    :type name: str

    """

    def __init__(
        self,
        fraction: float,
        mass_bins: Quantity["mass"],
        stars: Quantity["mass"],
        name="",
    ):
        self.fraction: float = fraction
        self.mass_bins: float = mass_bins.to(u.Msun).value
        self.stars: float = stars.to(u.Msun).value
        assert np.min(self.fraction) >= 0
        assert np.max(self.fraction) <= 1
        assert np.min(self.mass_bins) >= 0
        super().__init__(a=self.fraction, b=self.mass_bins, c=self.stars, name=name)


class StepFraction(Fraction):
    """
    A simple step function binary fraction, with a binary fraction
    of 0 below the changeover mass and 1 above the changeover mass

    :param mass: Changeover mass between binary fractions of 0 and 1
    :type mass: astropy mass unit

    """

    def __init__(
        self, fraction: float, mass_bins: Quantity["mass"], stars: Quantity["mass"]
    ):
        self.name = "Step"
        assert len(fraction) == 2
        super().__init__(fraction, mass_bins, stars, name=self.name)

    def binary_fraction(self) -> np.float64:
        """
        Binary fraction as a function of stellar mass
        for the step function binary fraction
        Returns the probability to be in a binary
        """
        prob = np.piecewise(
            self.stars,
            [self.stars < self.mass_bins, self.stars >= self.mass_bins],
            self.fraction,
        )
        return prob

    def sample(self) -> bool:
        """
        Determine which stars are primaries

        :return: Boolean array
        :rtype: bool
        """
        _sample = np.random.rand(len(self.stars))
        _binary = np.zeros(len(self.stars), dtype=bool)
        _select = np.where(_sample <= self.binary_fraction())
        _binary[_select] = np.ones(len(_select), dtype=bool)
        return _binary


class MassRatio(rv_continuous):
    """
    This class is the superclass of all mass ratio distributions
    TODO (CCC, 04/02/2025): q is currently independent of M1 and a

    :param min_q: Minimum mass ratio
    :type min_q: float
    :param max_q: Maximum mass ratio
    :type max_q: float
    :param stars: Primaries
    :type stars: Quantity["mass"], list of stellar masses
    :param name: Name for the scipy.stats rv_continuous instance
    :type name: str
    """

    def __init__(self, min_q: float, max_q: float, name=""):
        self.min_q = min_q
        self.max_q = max_q
        assert self.min_q > 0
        assert self.max_q >= min_q
        assert self.max_q <= 1
        super().__init__(a=self.min_q, b=self.max_q, name=name)

    def sample(self, N: int) -> Quantity["length"]:
        """
        :param N: Number of stars to draw
        :type N: int
        :return: List of semi-major axes of stars
        :rtype: Quantity["length"]
        """
        return self.rvs(size=N)


class UniformMassRatio(MassRatio):
    """
    A simple loguniform distribution of semi-major axes
    with lower and upper bounds
    """

    def __init__(self, min_q: float, max_q: float, name="uniform"):
        self.name = name
        super().__init__(min_q, max_q, name=self.name)

    def _pdf(self, x: np.float64, *args) -> np.float64:
        rv = uniform(self.min_q, self.max_q - self.min_q)
        # 2nd argument of scipy.stats.uniform is RANGE, not upper bound
        return rv.pdf(x)

    def _ppf(self, q: np.float64, *args) -> np.float64:
        rv = uniform(self.min_q, self.max_q - self.min_q)
        return rv.ppf(q)


class Period(rv_continuous):
    """
    This class is the superclass of all orbital period distributions

    :param min_p: Minimum orbital period
    :type min_p: astropy time unit
    :param max_p: Maximum orbital period
    :type max_p: astropy time unit
    :param stars: Primaries
    :type stars: Quantity["mass"], list of stellar masses
    :param name: Name for the scipy.stats rv_continuous instance
    :type name: str
    """

    def __init__(self, min_p: Quantity["time"], max_p: Quantity["time"], name=""):
        self.min_p: float = min_p.to(u.d).value
        self.max_p: float = max_p.to(u.d).value
        assert self.min_p > 0
        assert self.max_p > self.min_p
        super().__init__(a=self.min_p, b=self.max_p, name=name)

    def sample(self, N: int) -> Quantity["time"]:
        """
        :param N: Number of stars to draw
        :type N: int
        :return: List of orbital periods
        :rtype: Quantity["time"]
        """
        return self.rvs(size=N) * u.d


class LogUniformPeriod(Period):
    """
    A simple step function distribution of semi-major axes, with a log-uniform
    distribution above and below the changeover semi-major axis

    :param sma: Changeover semi-major axis
    :type sma: astropy length unit
    :param ratio: Ratio of probabilities for close and wide binaries
    :type ratio: float
    A simple loguniform distribution of semi-major axes

    """

    def __init__(
        self, min_p: Quantity["time"], max_p: Quantity["time"], name="loguniform"
    ):
        self.name = name
        super().__init__(min_p, max_p, name=self.name)

    def _pdf(self, x: np.float64, *args) -> np.float64:
        rv = loguniform(self.min_p, self.max_p)
        return rv.pdf(x)

    def _ppf(self, q: np.float64, *args) -> np.float64:
        rv = loguniform(self.min_p, self.max_p)
        return rv.ppf(q)


class Semimajor(rv_continuous):
    """
    This class is the superclass of all semi-major axis distributions

    :param min_a: Minimum semi-major axis
    :type min_a: astropy length unit
    :param max_a: Maximum semi-major axis
    :type max_a: astropy length unit
    :param stars: Primaries
    :type stars: Quantity["mass"], list of stellar masses
    :param name: Name for the scipy.stats rv_continuous instance
    :type name: str
    """

    def __init__(self, min_a: Quantity["length"], max_a: Quantity["length"], name=""):
        self.min_a: float = min_a.to(u.au).value
        self.max_a: float = max_a.to(u.au).value
        assert self.min_a > 0
        assert self.max_a > self.min_a
        super().__init__(a=self.min_a, b=self.max_a, name=name)

    def sample(self, N: int) -> Quantity["length"]:
        """
        :param N: Number of stars to draw
        :type N: int
        :return: List of semi-major axes of stars
        :rtype: Quantity["length"]
        """
        return self.rvs(size=N) * u.au


class LogUniformSemimajor(Semimajor):
    """
    A simple step function distribution of semi-major axes, with a log-uniform
    distribution above and below the changeover semi-major axis

    :param sma: Changeover semi-major axis
    :type sma: astropy length unit
    :param ratio: Ratio of probabilities for close and wide binaries
    :type ratio: float
    A simple loguniform distribution of semi-major axes
    """

    def __init__(
        self, min_a: Quantity["length"], max_a: Quantity["length"], name="loguniform"
    ):
        self.name = name
        super().__init__(min_a, max_a, name=self.name)

    def _pdf(self, x: np.float64, *args) -> np.float64:
        rv = loguniform(self.min_a, self.max_a)
        return rv.pdf(x)

    def _ppf(self, q: np.float64, *args) -> np.float64:
        rv = loguniform(self.min_a, self.max_a)
        return rv.ppf(q)
