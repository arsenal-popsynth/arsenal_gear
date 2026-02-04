# Formation context: the form_context variable accumulates
# everything needed to generate a stellar population.
imf = dist_funcs.imf.Salpeter(0.08 * u.Msun, 100 * u.Msun, alpha=2.3)
form_context = FormationContext(imf=imf, mass=1e6 * u.Msun, metals=0.01)

# Evolution context: the evolve_context variable accumulates
# everything needed to evolve a stellar population and get yields
explodability = lambda stars, t0, t1: np.logical_and(
    feedbacks.sn.explodable_mass_range(8 * u.Msun, 40 * u.Msun)(stars),
    feedbacks.sn.explodable_lifetime_range(t0, t1, feedbacks.sn.lifetimes_Raiteri)(
        stars
    ),
)

sn = feedbacks.SNFeedbackMechanism(
    energy_func=lambda stars: feedbacks.sn.constant_energy(stars, 1e51 * u.erg),
    mass_func=feedbacks.sn.massloss_Raiteri,
    metals_func=feedbacks.sn.metals_Raiteri,
    explodability_func=explodability,
)

evolve_context = EvolutionContext(mechanisms=[sn])

# out_times is a list of times at which to output results
out_times = np.logspace(5, 8, 10) * u.yr

# out_qtys is a list of quantities we want to output.
out_qtys = ["count", "energy", "mass", "metals_total", "Fe", "O"]
