imf = dist_funcs.imf.Salpeter(0.08 * u.Msun, 100 * u.Msun, alpha=2.3)

form_context = FormationContext(imf=imf, mass=1e6 * u.Msun, metals=0.01)

starpop = form_context.generate_population()

explodability = lambda stars, t0, t1: np.logical_and(
    feedbacks.sn.explodable_mass_range(8 * u.Msun, 40 * u.Msun)(stars),
    feedbacks.sn.explodable_lifetime_range(t0, t1, feedbacks.sn.lifetimes_Raiteri)(
        stars
    ),
)

sn_context = feedbacks.SNFeedbackMechanism(
    energy_func=lambda stars: feedbacks.sn.constant_energy(stars, 1e51 * u.erg),
    mass_func=feedbacks.sn.massloss_Raiteri,
    metals_func=feedbacks.sn.metals_Raiteri,
    explodability_func=explodability,
)

out_times = np.logspace(1, 8, 50) * u.yr

out_qtys = [
    lambda t0, t1: sn_context.count(starpop, t0, t1),
    lambda t0, t1: sn_context.energy(starpop, t0, t1),
    lambda t0, t1: sn_context.mass(starpop, t0, t1),
    lambda t0, t1: sn_context.metals_total(starpop, t0, t1),
    lambda t0, t1: sn_context.metals_species(starpop, "Fe", t0, t1),
    lambda t0, t1: sn_context.metals_species(starpop, "O", t0, t1),
]
