# Arsenal Gear Architecture Summary

## Overview

Four primary classes work together through composition:

- **`SynthPop`**: Main user-facing API that coordinates formation, evolution, and stellar products
- **`FormationContext`**: Encapsulates initial conditions and properties at birth
- **`EvolutionContext`**: Provides time-dependent stellar properties via `Isochrone` objects
- **`StellarProductsContext`**: Maps stellar properties to outputs (yields, mass-loss, spectra)

---

## Class Structure

### SynthPop

Coordinates the three contexts and handles population/time integration.

**Attributes**:
- `formation: FormationContext`
- `evolution: EvolutionContext`
- `products: StellarProductsContext`

**Public Methods**:
- Supernova: `nsn(t)`, `ndotsn(t)`
- Luminosity: `lbol(t)`, `lum(t, λ_min, λ_max)`, `mag(t, filter)`
- Per-star: `lbol_stars(t)`, `teff_stars(t)`, `mass_stars(t)`
- Feedback: `edot_wind(t)`, `edot_sn(t)`, `edot_total(t)`
- Mass loss: `mdot_wind(t)`, `mloss_wind(t)`, `mloss_sn(t)`, `mloss_total(t)`
- Yields: `yields(element, t)`, `yields_by_source(element, t)`

---

### FormationContext (Abstract Base)

Defines interface for atomic (single) populations.

**Concrete Implementations**:
- `SimpleFormationContext`: Single IMF, metallicity, age
- `BinaryFormationContext`: Includes binary parameters (future)
- `CompositeFormationContext`: Container for multiple atomic contexts (future)

**Abstract Interface**:
- `Mtot: Quantity["mass"]`
- `metallicity: float`
- `age: Quantity["time"]`
- `discrete: bool`
- `get_masses() -> Optional[Quantity["mass"]]`
- `get_mass_distribution() -> tuple`

---

### EvolutionContext

Provides stellar evolution data via `Isochrone` objects.

**Attributes**:
- `simple_evolution: IsochroneInterpolator`
- `binary_evolution: BinaryEvolutionModel` (future)

**Public Methods**:
- `get_isochrone(t) -> Isochrone`
- `get_lifetime(mass) -> Quantity["time"]`
- `get_mmax(t) -> Quantity["mass"]`
- `get_mmaxdot(t) -> Quantity["mass/time"]`

**Isochrone.qs contains arrays**:
- `initial_mass`, `current_mass`
- `log_Teff`, `log_L`, `log_g`
- `EEP`, `surface_abundances`, `mdot` (if available)

---

### StellarProductsContext

Maps stellar properties to outputs using configurable prescriptions.

**Attributes**:
- `source_config: Dict[str, OutputTypeConfig]`
- `yield_tables: Dict[str, Yields]`
- `mdot_prescriptions: Dict[str, MdotPrescription]`
- `atmosphere_models: Dict[str, Atmosphere]` (future)

**source_config Structure**:
```python
{
    'winds': {
        'phase_func': Callable[[Isochrone, str], np.ndarray],  # Returns bool mask
        'phases': {'MS': 'track', 'RGB': 'deJager1988', 'WR': 'NugisLamers2000'}
    },
    'yields': {
        'phase_func': Callable,
        'phases': {'CCSN': 'LC18', 'AGB': 'Karakas2010'}
    },
    'atmosphere': {
        'phase_func': Callable,
        'phases': {'hot': 'Kurucz', 'cool': 'MARCS'}
    }
}
```

**Public Methods** (all operate on `Isochrone`, return arrays):
- `get_phase_mask(iso, output_type, phase) -> np.ndarray`
- `get_mdot(iso) -> Quantity["mass/time"]`
- `get_edot_wind(iso) -> Quantity["power"]`
- `get_sn_yields(iso, element) -> Quantity["mass"]`
- `get_agb_yields(iso, element) -> Quantity["mass"]`
- `get_wind_yields(iso, element) -> Quantity["mass/time"]`
- `get_spectrum(iso) -> np.ndarray` (future)
- `get_ionizing_luminosity(iso) -> Quantity["photons/time"]` (future)

---

## Module Organization

```
arsenal_gear/
├── __init__.py                    # SynthPop
├── formation_context/
│   ├── __init__.py               # Abstract FormationContext
│   ├── simple.py                 # SimpleFormationContext
│   ├── binary.py                 # BinaryFormationContext (future)
│   └── composite.py              # CompositeFormationContext (future)
├── evolution_context.py           # EvolutionContext
├── stellar_products_context.py    # StellarProductsContext
├── dist_funcs/                    # IMF and binary distributions
├── stellar_evolution/             # IsochroneInterpolator, data structures
├── feedbacks/                     # Mass-loss prescriptions
├── element_yields/                # Yield tables
├── atmospheres/                   # Atmosphere models (future)
└── utils/
```

---

## Implementation Status

- **Phase 1** (Current): Simple Stellar Populations
- **Phase 2** (Future): Binary Populations
- **Phase 3** (Future): Composite Populations
- **Phase 4** (Future): Atmosphere Models
