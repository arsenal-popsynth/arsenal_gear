# Models Used in Existing Hydrocodes

## gasoline2/ChaNGa

Feedback from core collapse SN is governed based on a choice of IMF and
stellar lifetimes derived from Equation 3 of [Raiteri
1996](https://ui.adsabs.harvard.edu/abs/1996A%26A...315..105R/abstract).

These lifetimes depend on the mass *M* and metallicity *Z* of the star
(in solar units) and are given by:

$$t_*(\mathrm{yr}) = a_0(Z) + a_1(Z)\log(M) + a_2(Z)(\log(M))^2$$

With the metallicity-dependent coefficients *a*<sub>*i*</sub>(*Z*) given
by:

$$\begin{aligned}
a_0(Z) &= 10.13 + 0.07547\log(Z) - 0.008084\left(\log(Z)\right)^2 \\
a_1(Z) &= -4.424 - 0.7939\log(Z) - 0.1187\left(\log(Z)\right)^2 \\
a_2(Z) &= 1.262 + 0.3385\log(Z) + 0.05417\left(\log(Z)\right)^2
\end{aligned}$$
