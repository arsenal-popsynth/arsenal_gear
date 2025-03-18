==================================
Models Used in Existing Hydrocodes
==================================

gasoline2/ChaNGa
================
Feedback from core collapse SN is governed based on a choice of IMF and stellar 
lifetimes derived from Equation 3 of `Raiteri 1996 <https://ui.adsabs.harvard.edu/abs/1996A%26A...315..105R/abstract>`_.

These lifetimes depend on the mass :math:`M` and metallicity :math:`Z` of the star (in solar units) and are given by:

.. math::
    t_*(\mathrm{yr}) = a_0(Z) + a_1(Z)\log(M) + a_2(Z)\left(\log(M)\right)^2

With the metallicity-dependent coefficients :math:`a_i(Z)` given by:

.. math::
    a_0(Z) &= 10.13 + 0.07547\log(Z) - 0.008084\left(\log(Z)\right)^2 \\
    a_1(Z) &= -4.424 - 0.7939\log(Z) - 0.1187\left(\log(Z)\right)^2 \\
    a_2(Z) &= 1.262 + 0.3385\log(Z) + 0.05417\left(\log(Z)\right)^2