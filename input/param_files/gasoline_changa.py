# Formation context: the form_context variable accumulates
# everything needed to generate a stellar population.
imf = dist_funcs.imf.Salpeter(0.08 * u.Msun, 100 * u.Msun, alpha=2.3)
form_context = FormationContext(imf=imf, mass=1e6 * u.Msun, metals=0.01)


sn = feedbacks.SNFeedbackMechanism(
    energy_func=lambda stars: feedbacks.analytic.constant_energy(stars, 1e51 * u.erg),
    mass_func=feedbacks.analytic.mass_Raiteri,
    metals_func=feedbacks.analytic.metals_Raiteri,
    lifetime_func=feedbacks.analytic.lifetime_Raiteri,
    explodability_func=lambda stars: feedbacks.analytic.explodable_mass_range(
        8 * u.Msun, 40 * u.Msun
    ),
)
evolve_context = EvolutionContext(mechanisms=[sn])

# out_times is a list of times at which to output results
out_times = np.logspace(5, 8, 10) * u.yr

# out_qtys is a list of quantities we want to output.
out_qtys = ["count", "energy", "mass", "metals_total", "Fe", "O"]
