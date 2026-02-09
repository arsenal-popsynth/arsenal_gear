# Arsenal Gear Architecture

## Overview

Arsenal Gear's architecture is built around four primary classes that work together through composition:

- **`SynthPop`**: Main user-facing API that coordinates formation, evolution, and stellar products
- **`FormationContext`**: Encapsulates initial conditions and properties at birth
- **`EvolutionContext`**: Encapsulates time-dependent stellar evolution (stellar properties at time t)
- **`StellarProductsContext`**: Maps stellar properties to outputs (yields, mass-loss, spectra)

## Design Principles

### Composition Over Inheritance

The `SynthPop` class uses composition to integrate `FormationContext`, `EvolutionContext`, and `StellarProductsContext` instances, rather than inheriting from them. This provides:

- Flexibility to swap contexts independently
- Clear separation of concerns
- Easier testing of individual components
- Ability to reuse contexts across multiple populations

### Single Responsibility

- **`FormationContext`**: Knows about birth properties (mass, metallicity, IMF sampling)
- **`EvolutionContext`**: Knows about stellar evolution (stellar properties at time t via Isochrone)
- **`StellarProductsContext`**: Knows about mapping stellar properties to outputs (yields, feedback, spectra)
- **`SynthPop`**: Coordinates interactions and handles population-level integrations

## Class Structure

### SynthPop

**Purpose**: User-facing API that coordinates formation, evolution, and stellar products contexts.

**Key Responsibilities**:
- Create and manage `FormationContext`, `EvolutionContext`, and `StellarProductsContext` instances
- Provide polymorphic initialization (multiple ways to construct)
- Interface between the three contexts
- Handle population-level integrations (e.g., integrate stellar properties over IMF)
- Handle time integrations (cumulative yields, mass loss)
- Expose high-level methods for common queries

**Polymorphic Initialization**:

The `__init__()` method supports multiple construction patterns:

```python
# Option 1: Simple parameters (constructs contexts automatically)
pop = SynthPop(Mtot=1e6, metallicity=0.0, ...)

# Option 2: Explicit contexts (power users)
formation = FormationContext(...)
evolution = EvolutionContext(...)
products = StellarProductsContext(source_config={...})
pop = SynthPop(formation_context=formation, evolution_context=evolution,
               products_context=products)

# Option 3: Mix of both (provide some contexts, auto-construct others)
formation = FormationContext(...)
pop = SynthPop(formation_context=formation, metallicity=0.0, ...)
```

**Public Interface**:

The `SynthPop` class provides methods organized by physical process:

#### Supernova Methods
- `nsn(t)`: Cumulative number of supernovae by time t
  - Returns: Total number (or expected number) of SNe that have exploded by time t
  - Integrates over stars with mass > M_min_SN that have reached end of life
- `ndotsn(t)`: Supernova rate at time t (dN_SN/dt)
  - Returns: Instantaneous rate of supernovae at time t
  - Uses chain rule: dN/dt = (dN/dM) × (dM_max/dt)
  - This is fundamentally only defined for a continuous population,
    but we allow access even when you have a `discrete` `FormationContext`

#### Luminosity Methods - Population Integrated
- `lbol(t)`: Total bolometric luminosity at time t
  - Returns: Integrated luminosity over all stars in population
  - Discrete mode: sum over individual stars
  - Continuous mode: integrate L_bol(m,t) × IMF(m) dm
- `lum(t, wavelength_min, wavelength_max)`: Luminosity in wavelength range
  - Returns: Integrated luminosity between λ_min and λ_max
  - Requires spectral synthesis or filter convolution
- `mag(t, filter)`: Absolute magnitude in photometric filter
  - Returns: Population magnitude in specified filter (e.g., 'V', 'K', 'g', 'r')
  - Uses filter transmission curves and stellar spectra

#### Luminosity Methods - Per-Star
- `lbol_stars(t)`: Bolometric luminosity of each star at time t
  - Returns: Array of individual stellar luminosities
  - Only valid in discrete mode
  - Useful for finding individual luminous stars or performing custom integrations
- `lum_stars(t, wavelength_min, wavelength_max)`: Per-star luminosity in wavelength range
  - Returns: Array of luminosities for each star in specified wavelength range
  - Only valid in discrete mode
- `mag_stars(t, filter)`: Per-star absolute magnitude
  - Returns: Array of individual stellar magnitudes in specified filter
  - Only valid in discrete mode

#### Mechanical Feedback Methods - Population Integrated
- `edot_wind(t)`: Mechanical energy injection rate from stellar winds
  - Returns: dE/dt from mass loss in winds [erg/s or similar]
  - Includes contributions from main sequence, red giant, and Wolf-Rayet winds
- `edot_sn(t)`: Mechanical energy injection rate from supernovae
  - Returns: dE/dt from supernova explosions at time t
  - Typical SN energy: 10^51 erg (1 foe)
- `edot_total(t)`: Total mechanical energy injection rate
  - Returns: Combined energy rate from winds and supernovae
  - Equivalent to: edot_wind(t) + edot_sn(t)

#### Mechanical Feedback Methods - Per-Star
- `edot_wind_stars(t)`: Wind energy injection rate for each star
  - Returns: Array of dE/dt from winds for each individual star
  - Only valid in discrete mode
  - Zero for stars without significant winds
- `mdot_wind_stars(t)`: Wind mass loss rate for each star
  - Returns: Array of dM/dt from winds for each individual star
  - Only valid in discrete mode
  - Depends on stellar mass, evolutionary phase, and metallicity

#### Mass Loss Methods - Population Integrated
- `mdot_wind(t)`: Mass loss rate from stellar winds
  - Returns: dM/dt from winds at time t
  - Integrates mass loss rates over all living stars
- `mloss_wind(t)`: Cumulative mass lost to winds by time t
  - Returns: Total mass ejected via winds up to time t
  - Integral of mdot_wind(t') dt' from 0 to t
- `mloss_sn(t)`: Cumulative mass ejected in supernovae by time t
  - Returns: Total mass ejected via SN ejecta up to time t
  - Related to nsn(t) × M_ejecta (mass-dependent)
- `mloss_total(t)`: Total cumulative mass loss by time t
  - Returns: Combined mass loss from winds and supernovae

#### Chemical Yields Methods
- `yields(element, t)`: Cumulative elemental yields by time t
  - Returns: Total mass of specified element produced and ejected
  - Includes contributions from CCSN, SNIa, AGB, and winds
  - Example: `yields('O', t)` returns oxygen mass produced
- `yields_by_source(element, t)`: Yields broken down by source
  - Returns: Dictionary with keys ['ccsn', 'snia', 'agb', 'winds']
  - Useful for understanding dominant enrichment channels

#### Stellar Property Methods - Population Integrated
- `teff(t)`: Luminosity-weighted effective temperature at time t
  - Returns: <T_eff> weighted by L_bol for population
  - Useful for integrated spectral synthesis

#### Stellar Property Methods - Per-Star
- `teff_stars(t)`: Effective temperature of each star at time t
  - Returns: Array of T_eff for each individual star
  - Only valid in discrete mode
  - Useful for identifying hot stars, performing detailed spectral synthesis
- `mass_stars(t)`: Current mass of each star at time t
  - Returns: Array of current masses (accounts for mass loss)
  - Only valid in discrete mode
  - Different from initial masses due to winds and evolution

**Internal Workflow**:
1. Accept initialization parameters (polymorphic)
2. Construct `FormationContext` (atomic) if not provided
3. Construct `EvolutionContext` if not provided
4. For user queries:
   - Check if `formation` is atomic or composite
   - If atomic: compute directly using single metallicity/age
   - If composite: iterate over sub-populations and sum results
5. Perform population-level integrations using data from contexts

**Handling Atomic vs Composite Contexts**:
```python
def _dispatch(self, method_name: str, t: Quantity["time"], *args, **kwargs):
    """Generic dispatcher for atomic vs composite."""
    if isinstance(self.formation, CompositeFormationContext):
        # Sum over sub-populations
        results = []
        for subpop in self.formation.sub_populations:
            # Temporarily swap in atomic context
            original = self.formation
            self.formation = subpop
            result = getattr(self, f"_{method_name}_atomic")(t, *args, **kwargs)
            self.formation = original
            results.append(result)
        return sum(results)
    else:
        # Single atomic population
        return getattr(self, f"_{method_name}_atomic")(t, *args, **kwargs)
```

---

### FormationContext (Abstract Base Class for Atomic Populations)

**Purpose**: Defines the interface for **atomic** (single) population formation properties.

**Design Pattern**: Abstract base class for atomic contexts. `CompositeFormationContext` is separate (container pattern).

**Key Concept**: An atomic `FormationContext` represents a **single** population with:
- One metallicity value
- One formation time/age
- One IMF (for simple populations) or one set of binary parameters (for binary populations)

**Key Responsibilities**:
- Define common interface for atomic population properties
- Provide mass distributions for population-level integrations
- Support both discrete and continuous population modes
- Enable polymorphic use in `SynthPop`

**Abstract Interface** (all atomic subclasses must implement):

```python
class FormationContext(ABC):
    """Base class for atomic (single) stellar populations."""

    @property
    @abstractmethod
    def Mtot(self) -> Quantity["mass"]:
        """Total population mass"""

    @property
    @abstractmethod
    def metallicity(self) -> float:
        """Population metallicity (log10(Z/Zsun)) - single value for atomic populations"""

    @property
    @abstractmethod
    def age(self) -> Quantity["time"]:
        """Formation time / age of population - single value for atomic populations"""

    @property
    @abstractmethod
    def discrete(self) -> bool:
        """Whether this is a discrete or continuous population"""

    @abstractmethod
    def get_masses(self) -> Optional[Quantity["mass"]]:
        """Return stellar masses (discrete mode) or None (continuous mode)"""

    @abstractmethod
    def get_mass_distribution(self) -> tuple:
        """
        Return (masses, weights) for integration.

        For discrete mode: (sampled_masses, ones_array)
        For continuous mode: (mass_grid, imf_weights)
        """

    @abstractmethod
    def expected_stars_above(self, mass: Quantity["mass"]) -> float:
        """Expected number/fraction of stars above given mass"""
```

**Design Note**: This base class is for **atomic populations only**. The `CompositeFormationContext` does NOT inherit from this class - it's a separate container class.

---

### SimpleFormationContext

**Purpose**: Concrete implementation for Simple Stellar Populations (SSPs).

**Description**: Represents a single burst of star formation with one IMF, one metallicity, and one formation time. This is the initial implementation and most common use case.

**Attributes**:
- `Mtot`: Total population mass
- `metallicity`: Population metallicity (log10(Z/Zsun))
- `age`: Formation time (default: 0)
- `imf`: IMF object (e.g., Salpeter, Kroupa)
- `masses`: Sampled stellar masses (if discrete mode)
- `Nstar`: Number of stars (discrete) or expected number (continuous)
- `discrete`: Boolean flag for discrete vs. continuous mode

**Methods**:
- Implements all abstract methods from `FormationContext`
- `sample_masses()`: Sample stellar masses from IMF
- `resample()`: Re-sample masses with same IMF and Mtot

**Usage Example**:
```python
from arsenal_gear.dist_funcs.imf import Kroupa

formation = SimpleFormationContext(
    Mtot=1e6 * u.Msun,
    metallicity=0.0,
    age=0 * u.Myr,
    imf=Kroupa(0.08, 100),
    discrete=True
)
```

---

### BinaryFormationContext (Future)

**Purpose**: Concrete implementation for binary stellar populations.

**Description**: Represents a population including binary systems with period, mass ratio, and eccentricity distributions.

**Additional Attributes** (beyond SimpleFormationContext):
- `binary_fraction`: Fraction of stars in binaries
- `period_distribution`: Orbital period distribution
- `mass_ratio_distribution`: Mass ratio (q = M2/M1) distribution
- `eccentricity_distribution`: Orbital eccentricity distribution
- `primary_masses`: Primary star masses
- `secondary_masses`: Secondary star masses
- `periods`: Orbital periods
- `eccentricities`: Orbital eccentricities

**Additional Methods**:
- `sample_binaries()`: Sample binary parameters
- `get_primary_masses()`: Return primary star masses
- `get_secondary_masses()`: Return secondary star masses

**Note**: Binary populations will require extensions to `EvolutionContext` to handle binary evolution, mass transfer, and mergers.

---

### CompositeFormationContext (Future)

**Purpose**: Container for multiple atomic `FormationContext` instances.

**Design Pattern**: **Container/Composite pattern** - does NOT inherit from `FormationContext`.

**Description**: Combines multiple atomic `FormationContext` instances (simple or binary populations) to represent complex star formation histories with multiple bursts, metallicities, and/or formation times.

**Key Concept**: This is a **container**, not an atomic population. It holds multiple independent populations, each with their own metallicity and age.

**Attributes**:
- `sub_populations`: List of atomic `FormationContext` instances (Simple or Binary)
- `Mtot`: Total mass across all sub-populations (sum of component Mtots)
- `metallicities`: Array of metallicities (one per sub-population)
- `ages`: Array of formation times (one per sub-population)

**Methods**:
- `add_population(context: FormationContext)`: Add a sub-population
- `get_sub_population(index: int)`: Access individual sub-population
- `num_populations`: Number of sub-populations

**Usage Example**:
```python
# Multiple bursts at different times and metallicities
burst1 = SimpleFormationContext(Mtot=1e6*u.Msun, metallicity=0.0, age=0*u.Myr)
burst2 = SimpleFormationContext(Mtot=5e5*u.Msun, metallicity=-0.5, age=100*u.Myr)

composite = CompositeFormationContext(sub_populations=[burst1, burst2])

# Access properties
print(composite.Mtot)              # 1.5e6 Msun (sum)
print(composite.metallicities)     # [0.0, -0.5]
print(composite.ages)              # [0 Myr, 100 Myr]
print(composite.num_populations)   # 2
```

**How SynthPop Uses It**:
```python
# In SynthPop.lbol(t)
if isinstance(self.formation, CompositeFormationContext):
    # Iterate over sub-populations
    total = 0
    for subpop in self.formation.sub_populations:
        total += self._compute_lbol_atomic(subpop, t)
    return total
else:
    # Single atomic population
    return self._compute_lbol_atomic(self.formation, t)
```

**Use Cases**:
- Complex star formation histories (multiple bursts)
- Populations with metallicity spreads
- Globular clusters with multiple populations
- Galaxies with extended star formation
- Mixed simple + binary populations

---

### FormationContext Design Rationale

**Architecture Pattern: Atomic Contexts + Container**

This design uses two complementary patterns:

1. **Atomic Populations** (`FormationContext` base class):
   - `SimpleFormationContext` and `BinaryFormationContext` are atomic
   - Each has single metallicity, single age
   - Inherit from abstract `FormationContext` base class

2. **Composite Container** (`CompositeFormationContext`):
   - Does NOT inherit from `FormationContext`
   - Container holding multiple atomic contexts
   - Has array properties: `metallicities`, `ages`

**Why this approach?**

1. **Conceptual Clarity**: An atomic population *has* a metallicity (not an array of identical values)
2. **Type Correctness**: Atomic contexts return scalars; composite returns arrays
3. **Natural Composition**: Composite literally contains other contexts
4. **Clean Separation**: Each atomic context is fully self-contained
5. **Easy Querying**: "What's the metallicity?" has clear meaning for atomic populations
6. **Straightforward Iteration**: `SynthPop` can loop over sub-populations in composite case

**SynthPop Handles Both Cases**:

```python
def lbol(self, t):
    """Compute bolometric luminosity at time t."""
    if isinstance(self.formation, CompositeFormationContext):
        # Composite: sum over sub-populations
        return sum(
            self._lbol_atomic(subpop, t)
            for subpop in self.formation.sub_populations
        )
    else:
        # Atomic: single population
        return self._lbol_atomic(self.formation, t)

def _lbol_atomic(self, formation: FormationContext, t: Quantity["time"]):
    """Helper: compute Lbol for a single atomic population."""
    masses = formation.get_masses()
    metallicity = formation.metallicity  # Single value
    age = formation.age                  # Single value
    # ... compute luminosity using EvolutionContext ...
```

**Benefits of Container Pattern**:

1. **Simple populations stay simple**: No unnecessary array wrapping
2. **Type checking is clear**: `isinstance()` distinguishes atomic from composite
3. **Code duplication is minimal**: Atomic case is the base implementation
4. **Debugging is easier**: Single values vs arrays are immediately obvious
5. **Performance**: Direct access to scalar properties for common case

**Forward Compatibility**:

- Atomic `FormationContext` interface is stable
- New atomic types (e.g., `BinaryFormationContext`) follow same pattern
- `CompositeFormationContext` can hold any mix of atomic types
- `SynthPop` iteration logic handles all cases

**Implementation Strategy**:

1. **Phase 1** (Current): Implement `SimpleFormationContext` only
2. **Phase 2** (Future): Add `BinaryFormationContext` + binary evolution
3. **Phase 3** (Future): Add `CompositeFormationContext` container

**Data Flow**:
- Atomic contexts created directly or auto-constructed by `SynthPop`
- Composite contexts hold references to multiple atomic contexts
- `SynthPop` queries atomic contexts directly (simple) or iterates (composite)
- `EvolutionContext` always receives atomic context properties (single metallicity, single age)

---

### EvolutionContext

**Purpose**: Encapsulates stellar evolution models and provides time-dependent stellar properties.

**Design Pattern**: Single class with type dispatch - dispatches to appropriate evolution model based on `FormationContext` type.

**Key Responsibilities**:
- Dispatch to appropriate evolution model (simple or binary) based on formation context type
- Provide time-dependent stellar properties via `Isochrone` objects
- Provide stellar lifetimes and maximum living mass calculations

**Composed Components**:

```python
class EvolutionContext:
    # Evolution models (type-specific)
    simple_evolution: IsochroneInterpolator     # For SimpleFormationContext
    binary_evolution: BinaryEvolutionModel      # For BinaryFormationContext (future)
```

**Public Interface**:

```python
class EvolutionContext:
    def __init__(self, metallicity: float, ...):
        """Initialize with metallicity and load appropriate stellar models."""

    def get_isochrone(self, t: Quantity["time"]) -> Isochrone:
        """
        Returns Isochrone with arrays of stellar properties at time t.

        The Isochrone.qs dict contains arrays:
        - initial_mass, current_mass
        - log_Teff, log_L, log_g
        - EEP (Equivalent Evolutionary Phase)
        - surface_abundances (if available)
        - mdot (track-based mass-loss, if available)
        """

    def get_lifetime(self, mass: Quantity["mass"]) -> Quantity["time"]:
        """Returns lifetime for a star of given initial mass."""

    def get_mmax(self, t: Quantity["time"]) -> Quantity["mass"]:
        """Returns maximum initial mass still alive at time t."""

    def get_mmaxdot(self, t: Quantity["time"]) -> Quantity["mass/time"]:
        """Returns rate of change of maximum living mass."""
```

**Isochrone Structure** (from `stellar_evolution/se_data_structures.py`):

The `Isochrone` dataclass stores arrays of stellar properties, enabling vectorized operations:

```python
@dataclass
class Isochrone:
    age: Quantity["time"]
    eep_name: str       # Column name for EEP
    mini_name: str      # Column name for initial mass
    lteff_name: str     # Column name for log(Teff)
    llbol_name: str     # Column name for log(Lbol)
    qs: dict            # Dictionary of property arrays

    # qs contains arrays like:
    # - qs['initial_mass']: array of initial masses
    # - qs['log_Teff']: array of effective temperatures
    # - qs['log_L']: array of bolometric luminosities
    # - qs['EEP']: array of evolutionary phases
    # - qs['log_g']: surface gravity (if available)
    # - qs['mdot']: track-based mass-loss rates (if available)
```

**Design Rationale**:

1. **Returns Isochrone objects**: Enables vectorized operations downstream
2. **Type dispatch**: `EvolutionContext` internally routes to appropriate evolution model
3. **Separation from outputs**: Yields, feedback, and atmospheres handled by `StellarProductsContext`
4. **Track data included**: If tracks provide mass-loss rates or surface abundances, these are included in the Isochrone

**Data Dependencies**:
- Metallicity set at initialization
- Evolution models loaded based on metallicity
- Binary evolution additionally depends on orbital parameters (future)

---

### StellarProductsContext

**Purpose**: Maps stellar properties to physical outputs (yields, mass-loss rates, spectra) using configurable prescriptions.

**Design Pattern**: Configurable routing via `source_config` with phase-based mask functions for vectorized operations.

**Key Responsibilities**:
- Apply yield prescriptions (CCSN, AGB, wind yields)
- Apply mass-loss rate prescriptions (MS, RGB, WR winds)
- Apply atmosphere models for spectral synthesis (future)
- Route to track-based values or external prescriptions based on configuration

**Composed Components**:

```python
class StellarProductsContext:
    # Configuration for routing to sources
    source_config: Dict[str, OutputTypeConfig]

    # Prescription instances (initialized based on source_config)
    yield_tables: Dict[str, Yields]
    mdot_prescriptions: Dict[str, MdotPrescription]
    atmosphere_models: Dict[str, Atmosphere]  # Future
```

**source_config Structure**:

The `source_config` dictionary specifies which prescription to use for each output type and evolutionary phase:

```python
source_config = {
    'winds': {
        'phase_func': classify_wind_phase,  # (Isochrone, phase: str) -> np.ndarray[bool]
        'phases': {
            'MS': 'track',              # Use track-based mass-loss
            'RGB': 'deJager1988',        # Use de Jager prescription
            'WR': 'NugisLamers2000',     # Use Nugis & Lamers prescription
        }
    },
    'yields': {
        'phase_func': classify_yield_phase,
        'phases': {
            'CCSN': 'LC18',              # Limongi & Chieffi 2018
            'AGB': 'Karakas2010',        # Karakas 2010
        }
    },
    'atmosphere': {
        'phase_func': classify_atm_phase,
        'phases': {
            'hot': 'Kurucz',
            'cool': 'MARCS',
        }
    }
}
```

**Phase Functions Return Masks**:

Phase functions take an `Isochrone` and phase name, returning a boolean mask for vectorized operations:

```python
def classify_wind_phase(iso: Isochrone, phase: str) -> np.ndarray:
    """Returns boolean mask for stars in given wind phase."""
    EEP = iso.qs[iso.eep_name]
    log_Teff = iso.qs[iso.lteff_name]
    log_L = iso.qs[iso.llbol_name]

    if phase == 'MS':
        return (EEP >= 202) & (EEP < 454)
    elif phase == 'RGB':
        return (EEP >= 454) & (EEP < 605) & (log_Teff < 3.7)
    elif phase == 'WR':
        return (log_Teff > 4.5) & (log_L > 5.0)
    else:
        return np.zeros(len(EEP), dtype=bool)
```

**Public Interface**:

All methods operate on `Isochrone` objects and return arrays:

```python
class StellarProductsContext:
    def __init__(self, source_config: dict):
        """Initialize and create prescription instances from config."""

    # --- Phase Masks ---
    def get_phase_mask(self, iso: Isochrone, output_type: str, phase: str) -> np.ndarray:
        """Returns boolean mask for stars in given phase."""
        phase_func = self.source_config[output_type]['phase_func']
        return phase_func(iso, phase)

    # --- Mass-Loss Rates (Vectorized) ---
    def get_mdot(self, iso: Isochrone) -> Quantity["mass/time"]:
        """
        Returns array of mass-loss rates for all stars in isochrone.
        Routes to track or prescription based on phase config.
        """
        mdot = np.zeros(len(iso.qs[iso.mini_name])) * u.Msun/u.yr

        for phase, source in self.source_config['winds']['phases'].items():
            mask = self.get_phase_mask(iso, 'winds', phase)
            if source == 'track':
                mdot[mask] = iso.qs['mdot'][mask]
            else:
                mdot[mask] = self.mdot_prescriptions[source].compute(iso, mask)

        return mdot

    def get_edot_wind(self, iso: Isochrone) -> Quantity["power"]:
        """Returns array of wind mechanical energy rates."""

    # --- Yields (Vectorized) ---
    def get_sn_yields(self, iso: Isochrone, element: str) -> Quantity["mass"]:
        """
        Returns array of SN yields for given element.
        Non-zero only for stars that will explode as SNe.
        """

    def get_agb_yields(self, iso: Isochrone, element: str) -> Quantity["mass"]:
        """Returns array of AGB yields for given element."""

    def get_wind_yields(self, iso: Isochrone, element: str) -> Quantity["mass/time"]:
        """Returns array of wind yield rates for given element."""

    # --- Atmosphere/Spectra (Future) ---
    def get_spectrum(self, iso: Isochrone) -> np.ndarray:
        """Returns spectra for all stars in isochrone."""

    def get_ionizing_luminosity(self, iso: Isochrone) -> Quantity["photons/time"]:
        """Returns ionizing photon rates for all stars."""
```

**Internal Routing Logic**:

```python
def get_mdot(self, iso: Isochrone) -> Quantity["mass/time"]:
    """Vectorized mass-loss rate computation with phase routing."""
    n_stars = len(iso.qs[iso.mini_name])
    mdot = np.zeros(n_stars) * u.Msun/u.yr

    for phase, source in self.source_config['winds']['phases'].items():
        # Get mask for this phase
        mask = self.get_phase_mask(iso, 'winds', phase)

        if not np.any(mask):
            continue

        # Route to appropriate source
        if source == 'track':
            # Use track-based mass-loss from isochrone
            mdot[mask] = iso.qs['mdot'][mask]
        else:
            # Use external prescription
            mdot[mask] = self.mdot_prescriptions[source].compute(iso, mask)

    return mdot
```

**Design Rationale**:

1. **Vectorized operations**: All methods operate on arrays via `Isochrone`, not individual stars
2. **Mask-based phase selection**: Efficient boolean indexing for phase-specific calculations
3. **Configurable routing**: `source_config` allows flexible choice between track data and prescriptions
4. **`'track'` as special source**: Indicates data should come from the `Isochrone.qs` dict
5. **Separation from evolution**: Cleanly separates "what is the stellar state" from "what does the star produce"

**Prompt vs. Continuous Feedback**:

The `StellarProductsContext` provides **instantaneous** quantities:
- **Continuous** (winds): `get_mdot()`, `get_edot_wind()`, `get_wind_yields()` - rates at time t
- **Prompt** (SNe, AGB): `get_sn_yields()`, `get_agb_yields()` - yields from SN up to this time

Time integration (cumulative yields, total mass loss) is handled by `SynthPop`.

**Module Integration**:

| Module | Location | Purpose |
|--------|----------|---------|
| `Yields` implementations | `element_yields/` | CCSN, AGB, wind yield tables |
| `MdotPrescription` implementations | `feedbacks/` (new) | Mass-loss rate prescriptions |
| `Atmosphere` implementations | `atmospheres/` (future) | Spectral synthesis |

---

## Data Flow

### Initialization Flow

#### For Atomic Population:
```
User Input → SynthPop.__init__()
             ├→ Create/receive atomic FormationContext
             │  ├→ Initialize IMF
             │  └→ Sample masses (if discrete)
             ├→ Create/receive EvolutionContext
             │  └→ Initialize IsochroneInterpolator (or BinaryEvolutionModel)
             └→ Create/receive StellarProductsContext
                ├→ Parse source_config
                ├→ Initialize yield tables (e.g., LimongiChieffi2018)
                ├→ Initialize mdot prescriptions (e.g., VinkMdot)
                └→ Initialize atmosphere models (future)
```

#### For Composite Population:
```
User Input → SynthPop.__init__()
             ├→ Create/receive CompositeFormationContext
             │  └→ Contains multiple atomic FormationContext instances
             │     ├→ Each has its own IMF, metallicity, age
             │     └→ Each has sampled masses (if discrete)
             ├→ Create/receive EvolutionContext
             │  └→ Initialize stellar evolution models
             └→ Create/receive StellarProductsContext
                └→ Shared across all sub-populations

Note: EvolutionContext and StellarProductsContext are shared across all sub-populations.
```

### Query Flow (Example: `SynthPop.lbol(t)`)

#### For Atomic FormationContext (Simple or Binary):
```
User calls: pop.lbol(t)
    ↓
SynthPop.lbol(t):
    ├→ Check: isinstance(formation, CompositeFormationContext)? No
    ├→ Get isochrone: iso = evolution.get_isochrone(t)
    │   ├→ If Simple: IsochroneInterpolator returns Isochrone
    │   └→ If Binary: BinaryEvolutionModel returns Isochrone (accounts for mass transfer, mergers)
    ├→ Extract L_bol array: lbol_array = 10**iso.qs['log_L']
    ├→ If discrete mode:
    │   └→ Interpolate to sampled masses and sum
    └→ If continuous mode:
        └→ Integrate over IMF: ∫ L_bol(m) × IMF(m) dm
```

#### For CompositeFormationContext:
```
User calls: pop.lbol(t)
    ↓
SynthPop.lbol(t):
    ├→ Check: isinstance(formation, CompositeFormationContext)? Yes
    ├→ Initialize: total_lbol = 0
    ├→ For each subpop in formation.sub_populations:
    │   ├→ Get isochrone for subpop's age relative to t
    │   │   ├→ If Simple: IsochroneInterpolator returns Isochrone
    │   │   └→ If Binary: BinaryEvolutionModel returns Isochrone
    │   ├→ Compute subpop luminosity (discrete sum or IMF integral)
    │   └→ Add to total_lbol
    └→ Return total_lbol

Note: EvolutionContext handles mixed populations
      (e.g., burst1 = SimpleFormationContext, burst2 = BinaryFormationContext)
```

### Query Flow (Example: `SynthPop.mdot_wind(t)`)

```
User calls: pop.mdot_wind(t)
    ↓
SynthPop.mdot_wind(t):
    ├→ Get isochrone: iso = evolution.get_isochrone(t)
    ├→ Get mass-loss rates: mdot_array = products.get_mdot(iso)
    │   ↓
    │   StellarProductsContext.get_mdot(iso):
    │       ├→ For each phase in source_config['winds']['phases']:
    │       │   ├→ Get mask: mask = phase_func(iso, phase)
    │       │   ├→ If source == 'track':
    │       │   │   └→ mdot[mask] = iso.qs['mdot'][mask]
    │       │   └→ Else:
    │       │       └→ mdot[mask] = prescription.compute(iso, mask)
    │       └→ Return mdot array
    │   ↓
    ├→ If discrete mode:
    │   └→ Interpolate to sampled masses and sum
    └→ If continuous mode:
        └→ Integrate over IMF: ∫ ṁ(m) × IMF(m) dm
```

### Query Flow (Example: `SynthPop.yields('O', t)`)

```
User calls: pop.yields('O', t)
    ↓
SynthPop.yields('O', t):
    ├→ Compute cumulative SN yields (prompt):
    │   ├→ Find stars that have died by time t (mass > mmax(t))
    │   ├→ For dead stars: yields_sn = products.get_sn_yields(iso, 'O')
    │   └→ Sum over dead stars
    │
    ├→ Compute cumulative wind yields (continuous):
    │   ├→ For each time step t' from 0 to t:
    │   │   ├→ iso = evolution.get_isochrone(t')
    │   │   └→ wind_yield_rate = products.get_wind_yields(iso, 'O')
    │   └→ Integrate: ∫₀ᵗ wind_yield_rate(t') dt'
    │
    └→ Return total = sn_yields + wind_yields
```

### Integration Pattern

For population-level properties, `SynthPop` handles integration:

1. **Discrete Mode**: Interpolate isochrone to sampled masses and sum
   ```python
   # Get isochrone arrays
   iso = evolution.get_isochrone(t)
   lbol_grid = 10**iso.qs['log_L']
   mass_grid = iso.qs['initial_mass']

   # Interpolate to sampled masses
   lbol_sampled = np.interp(formation.masses, mass_grid, lbol_grid)
   total_lbol = np.sum(lbol_sampled)
   ```

2. **Continuous Mode**: Integrate over IMF
   ```python
   iso = evolution.get_isochrone(t)
   lbol_grid = 10**iso.qs['log_L']
   mass_grid = iso.qs['initial_mass']

   # Integrate L_bol(m) × IMF(m) dm
   total_lbol = integrate(lbol_grid * formation.imf.pdf(mass_grid), mass_grid)
   ```

---

## Module Organization

### Current Module Structure

```
arsenal_gear/
├── __init__.py                    # SynthPop (main API)
├── formation_context/             # Formation context classes (NEW)
│   ├── __init__.py               # Abstract FormationContext base
│   ├── simple.py                 # SimpleFormationContext
│   ├── binary.py                 # BinaryFormationContext (future)
│   └── composite.py              # CompositeFormationContext (future)
├── evolution_context.py           # EvolutionContext class (NEW)
├── stellar_products_context.py    # StellarProductsContext class (NEW)
├── population.py                  # StarPopulation, BinaryPopulation (data containers)
├── dist_funcs/                    # IMF and binary distributions
│   ├── imf.py                    # Used by FormationContext
│   └── binaries.py
├── stellar_evolution/            # Used by EvolutionContext
│   ├── isochrone.py
│   ├── se_data_structures.py
│   ├── data_reader.py
│   └── popsynth.py
├── feedbacks/                    # Used by StellarProductsContext
│   ├── __init__.py
│   ├── sn.py
│   └── mdot_prescriptions.py     # Mass-loss rate prescriptions (NEW)
├── element_yields/               # Used by StellarProductsContext
│   ├── yields.py
│   ├── yieldtables.py
│   └── [various implementations]
├── atmospheres/                  # Used by StellarProductsContext (future)
│   └── [atmosphere model implementations]
└── utils/                        # Shared utilities
```

### Dependency Graph

```
SynthPop
├── FormationContext (one of):
│   ├── Atomic FormationContext (abstract base):
│   │   ├── SimpleFormationContext
│   │   │   └── IMF (from dist_funcs/)
│   │   └── BinaryFormationContext (future)
│   │       ├── IMF (from dist_funcs/)
│   │       └── BinaryDistribution (from dist_funcs/)
│   │
│   └── CompositeFormationContext (container, future)
│       └── List[FormationContext] (atomic contexts only)
│
├── EvolutionContext
│   └── Evolution Models:
│       ├── simple_evolution: IsochroneInterpolator (from stellar_evolution/)
│       └── binary_evolution: BinaryEvolutionModel (from stellar_evolution/) - future
│
└── StellarProductsContext
    ├── source_config: Dict[str, OutputTypeConfig]
    ├── yield_tables: Dict[str, Yields] (from element_yields/)
    ├── mdot_prescriptions: Dict[str, MdotPrescription] (from feedbacks/)
    └── atmosphere_models: Dict[str, Atmosphere] (from atmospheres/) - future

Data Flow:
  FormationContext → provides masses, metallicity
  EvolutionContext → provides Isochrone (stellar properties at time t)
  StellarProductsContext → maps Isochrone to outputs (yields, mass-loss, spectra)
  SynthPop → integrates over population and time

Note: CompositeFormationContext does NOT inherit from FormationContext.
      It's a separate container class that holds atomic FormationContext instances.
```

---

## Design Decisions & Rationale

### Why Composition?

- **Modularity**: Contexts can be developed and tested independently
- **Reusability**: Same formation context can be used with different evolution models
- **Extensibility**: New context types (e.g., `BinaryFormationContext`) can be added without modifying core classes

### Why Separate EvolutionContext and StellarProductsContext?

The separation reflects a fundamental distinction in the physics:

1. **EvolutionContext**: "What is the stellar state at time t?"
   - Time-dependent stellar structure (L_bol, T_eff, log_g, mass, etc.)
   - Comes from stellar evolution models (isochrones, tracks)
   - Returns `Isochrone` objects with arrays of stellar properties

2. **StellarProductsContext**: "Given those properties, what does the star produce?"
   - Maps stellar properties to outputs (yields, mass-loss, spectra)
   - Configurable prescriptions per evolutionary phase
   - Can use track data OR external prescriptions

**Benefits of this separation**:
- **Swappable components**: Use different yield tables with same stellar evolution
- **Configurable routing**: Choose track-based or prescription-based outputs per phase
- **Vectorized operations**: All methods operate on `Isochrone` arrays
- **Clear data flow**: Evolution → Properties → Products

### Why Keep SynthPop as Coordinator?

- **User Convenience**: Simple, high-level API for common use cases
- **Backward Compatibility**: Existing code using `SynthPop` can continue to work
- **Integration Logic**: Population-level calculations (sums, integrals) naturally belong at this level
- **Time Integration**: Handles prompt (SN) vs. continuous (winds) feedback integration
- **Context Isolation**: Keeps `FormationContext`, `EvolutionContext`, and `StellarProductsContext` decoupled

### Polymorphic Initialization Benefits

- **Ease of Use**: Beginners can use simple parameters
- **Power User Control**: Advanced users can configure contexts explicitly
- **Gradual Learning Curve**: Users can start simple and grow into complexity

---

## Future Extensions

### Planned Enhancements

1. **BinaryFormationContext**: Extend formation context for binary populations
2. **Custom Evolution Models**: Pluggable evolution contexts for different stellar models
3. **Time-Dependent Metallicity**: Allow metallicity to vary with time (chemical evolution)
4. **Parallel Populations**: Manage multiple populations with different contexts
5. **Caching**: Cache expensive calculations in contexts

### Extensibility Points

- Add new `Yields` implementations for different yield tables (in `element_yields/`)
- Add new `MdotPrescription` implementations for different mass-loss prescriptions (in `feedbacks/`)
- Add new `Atmosphere` implementations for different spectral models (in `atmospheres/`)
- Add new `IMF` classes for different mass functions (in `dist_funcs/`)
- Add new phase classification functions for different evolutionary phase schemes
- Create custom `EvolutionContext` subclasses for different stellar evolution codes

---

## Implementation Status

### Phase 1: Simple Stellar Populations (Current Priority)
- [ ] Create abstract `FormationContext` base class
- [ ] Create `SimpleFormationContext` implementation
- [ ] Create `EvolutionContext` class (stellar evolution only)
- [ ] Create `StellarProductsContext` class with `source_config` structure
- [ ] Implement phase mask functions for wind/yield/atmosphere routing
- [ ] Refactor `SynthPop` to use composition with all three contexts
- [ ] Implement polymorphic initialization
- [ ] Update tests for new architecture
- [ ] Update documentation and examples
- [ ] Ensure backward compatibility

### Phase 2: Binary Populations (Future)
- [ ] Create `BinaryFormationContext` class
- [ ] Add `BinaryEvolutionModel` to `EvolutionContext`
- [ ] Ensure `StellarProductsContext` works with binary isochrones
- [ ] Add binary-specific methods to `SynthPop`
- [ ] Integrate binary evolution tracks
- [ ] Add tests for binary populations

### Phase 3: Composite Populations (Future)
- [ ] Create `CompositeFormationContext` class
- [ ] Handle multiple metallicities and ages in `SynthPop`
- [ ] Support complex star formation histories
- [ ] Add tests for composite populations

### Phase 4: Atmosphere Models (Future)
- [ ] Create `Atmosphere` base class
- [ ] Implement Kurucz atmosphere model
- [ ] Implement MARCS atmosphere model
- [ ] Add spectral synthesis methods to `StellarProductsContext`
- [ ] Add `get_spectrum()`, `get_ionizing_luminosity()` to `SynthPop`

---

## Notes

This architecture document is a living document and will be updated as the implementation progresses. See commit history for evolution of design decisions.
