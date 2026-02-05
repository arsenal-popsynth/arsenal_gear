# Arsenal Gear Architecture

## Overview

Arsenal Gear's architecture is built around three primary classes that work together through composition:

- **`SynthPop`**: Main user-facing API that coordinates formation and evolution
- **`FormationContext`**: Encapsulates initial conditions and properties at birth
- **`EvolutionContext`**: Encapsulates time-dependent evolution and physical processes

## Design Principles

### Composition Over Inheritance

The `SynthPop` class uses composition to integrate `FormationContext` and `EvolutionContext` instances, rather than inheriting from them. This provides:

- Flexibility to swap contexts independently
- Clear separation of concerns
- Easier testing of individual components
- Ability to reuse contexts across multiple populations

### Single Responsibility

- **`FormationContext`**: Knows about birth properties (mass, metallicity, IMF sampling)
- **`EvolutionContext`**: Knows about physical processes (stellar evolution, yields, feedback)
- **`SynthPop`**: Coordinates interactions and handles population-level integrations

## Class Structure

### SynthPop

**Purpose**: User-facing API that coordinates formation and evolution contexts.

**Key Responsibilities**:
- Create and manage `FormationContext` and `EvolutionContext` instances
- Provide polymorphic initialization (multiple ways to construct)
- Interface between the two contexts
- Handle population-level integrations (e.g., integrate stellar properties over IMF)
- Expose high-level methods for common queries

**Polymorphic Initialization**:

The `__init__()` method supports multiple construction patterns:

```python
# Option 1: Simple parameters (constructs contexts automatically)
pop = SynthPop(Mtot=1e6, metallicity=0.0, ...)

# Option 2: Explicit contexts (power users)
formation = FormationContext(...)
evolution = EvolutionContext(...)
pop = SynthPop(formation_context=formation, evolution_context=evolution)

# Option 3: Mix of both (provide one context, auto-construct the other)
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

**Purpose**: Encapsulates all time-dependent physical processes and evolution models for both simple and binary populations.

**Design Pattern**: Single class with type dispatch - dispatches to appropriate evolution model based on `FormationContext` type.

**Key Responsibilities**:
- Dispatch to appropriate evolution model (simple or binary) based on formation context type
- Provide time-dependent stellar properties (isochrone interpolation for simple, binary evolution for binaries)
- Handle feedback mechanisms (supernovae, winds, radiation) for both population types
- Compute elemental yields for both population types

**Composed Components**:

The `EvolutionContext` contains instances of specialized modules:

```python
class EvolutionContext:
    # Evolution models (type-specific)
    simple_evolution: IsochroneInterpolator     # For SimpleFormationContext
    binary_evolution: BinaryEvolutionModel      # For BinaryFormationContext (future)

    # Shared modules (work for both simple and binary)
    yields: Yields                              # Chemical enrichment
    feedback: FBMechanism                       # Energy and momentum feedback
```

**Module Integration**:

| Module | Location | Population Type | Purpose |
|--------|----------|-----------------|---------|
| `IsochroneInterpolator` | `stellar_evolution/` | Simple | Single star evolution (L_bol, T_eff, lifetimes) |
| `BinaryEvolutionModel` | `stellar_evolution/` (future) | Binary | Binary evolution (mass transfer, mergers, etc.) |
| `Yields` | `element_yields/` | Both | Chemical enrichment (CCSN, SNIa, AGB, winds) |
| `FBMechanism` | `feedbacks/` | Both | Energy and momentum feedback |

**Public Interface**:

Methods accept a `FormationContext` and dispatch internally:

```python
def lbol(self, formation: FormationContext, t: Quantity["time"]) -> Quantity["power"]:
    """Compute bolometric luminosity for given formation context at time t."""
    if isinstance(formation, SimpleFormationContext):
        return self._lbol_simple(formation, t)
    elif isinstance(formation, BinaryFormationContext):
        return self._lbol_binary(formation, t)
    else:
        raise TypeError(f"Unknown formation context type: {type(formation)}")

def teff(self, formation: FormationContext, t: Quantity["time"]) -> Quantity["temperature"]:
    """Compute effective temperature for given formation context at time t."""
    # Similar dispatch pattern

def yields_at_time(self, formation: FormationContext, element: str, t: Quantity["time"]) -> Quantity["mass"]:
    """Compute cumulative elemental yields for given formation context."""
    # Uses yields module (shared for both types)
```

**Internal Implementation** (Simple Populations):

```python
def _lbol_simple(self, formation: SimpleFormationContext, t: Quantity["time"]):
    """Compute Lbol for simple stellar population."""
    masses = formation.get_masses()
    metallicity = formation.metallicity
    age = formation.age

    # Use IsochroneInterpolator
    stellar_lbol = self.simple_evolution.lbol(masses, metallicity, t)
    return stellar_lbol

def _teff_simple(self, formation: SimpleFormationContext, t: Quantity["time"]):
    """Compute Teff for simple stellar population."""
    masses = formation.get_masses()
    metallicity = formation.metallicity

    stellar_teff = self.simple_evolution.teff(masses, metallicity, t)
    return stellar_teff
```

**Internal Implementation** (Binary Populations - Future):

```python
def _lbol_binary(self, formation: BinaryFormationContext, t: Quantity["time"]):
    """Compute Lbol for binary stellar population."""
    primaries = formation.get_primary_masses()
    secondaries = formation.get_secondary_masses()
    periods = formation.periods
    metallicity = formation.metallicity
    age = formation.age

    # Use BinaryEvolutionModel
    # Accounts for mass transfer, mergers, etc.
    binary_lbol = self.binary_evolution.lbol(primaries, secondaries, periods, metallicity, t)
    return binary_lbol
```

**Yields and Feedback** (Shared Logic):

The yields and feedback modules work for both simple and binary populations:

```python
def compute_yields(self, formation: FormationContext, element: str, t: Quantity["time"]):
    """Compute yields - works for both simple and binary."""
    # Get appropriate masses based on formation type
    if isinstance(formation, SimpleFormationContext):
        masses = formation.get_masses()
    elif isinstance(formation, BinaryFormationContext):
        # For binaries, consider both components
        masses = self._get_all_binary_masses(formation)

    # Use shared yields module
    return self.yields.compute(masses, formation.metallicity, element, t)
```

**Design Rationale**:

1. **Single EvolutionContext**: Users create one context that handles all population types
2. **Type dispatch**: `EvolutionContext` internally routes to appropriate evolution model
3. **Shared modules**: Yields and feedback logic reused for both simple and binary
4. **Matches FormationContext pattern**: Consistent with how `SynthPop` dispatches on atomic vs composite
5. **Works with CompositeFormationContext**: Since composite contains mixed simple/binary populations, a single `EvolutionContext` handles all

**Data Dependencies**:
- Receives `FormationContext` (atomic only) from `SynthPop`
- Extracts masses, metallicity, age from formation context
- Evolution models are metallicity-dependent
- Binary evolution additionally depends on orbital parameters

---

## Data Flow

### Initialization Flow

#### For Atomic Population:
```
User Input → SynthPop.__init__()
             ├→ Create/receive atomic FormationContext
             │  ├→ Initialize IMF
             │  └→ Sample masses (if discrete)
             └→ Create/receive EvolutionContext
                ├→ Initialize IsochroneInterpolator (simple_evolution)
                ├→ Initialize BinaryEvolutionModel (binary_evolution) - future
                ├→ Initialize Yields
                └→ Initialize FBMechanism
```

#### For Composite Population:
```
User Input → SynthPop.__init__()
             ├→ Create/receive CompositeFormationContext
             │  └→ Contains multiple atomic FormationContext instances
             │     ├→ Each has its own IMF, metallicity, age
             │     └→ Each has sampled masses (if discrete)
             └→ Create/receive EvolutionContext
                ├→ Initialize IsochroneInterpolator (simple_evolution)
                ├→ Initialize BinaryEvolutionModel (binary_evolution) - future
                ├→ Initialize Yields (shared)
                └→ Initialize FBMechanism (shared)

Note: EvolutionContext is shared across all sub-populations.
      It dispatches to appropriate evolution model based on each sub-population's type.
      Yields and feedback are shared modules that work for both simple and binary.
```

### Query Flow (Example: `SynthPop.lbol(t)`)

#### For Atomic FormationContext (Simple or Binary):
```
User calls: pop.lbol(t)
    ↓
SynthPop.lbol(t):
    ├→ Check: isinstance(formation, CompositeFormationContext)? No
    ├→ Call: _lbol_atomic(formation, t)
    │   ├→ Call EvolutionContext.lbol(formation, t)
    │   │   ├→ Check formation type: SimpleFormationContext or BinaryFormationContext?
    │   │   ├→ If Simple:
    │   │   │   ├→ Get masses from formation.get_masses()
    │   │   │   ├→ Get metallicity, age from formation
    │   │   │   └→ Call simple_evolution.lbol(masses, metallicity, t)
    │   │   │       └→ IsochroneInterpolator.lbol(masses, t)
    │   │   └→ If Binary:
    │   │       ├→ Get primary/secondary masses, periods from formation
    │   │       ├→ Get metallicity, age from formation
    │   │       └→ Call binary_evolution.lbol(primaries, secondaries, periods, metallicity, t)
    │   └→ Integrate over mass distribution (if continuous mode)
    └→ Return total luminosity
```

#### For CompositeFormationContext:
```
User calls: pop.lbol(t)
    ↓
SynthPop.lbol(t):
    ├→ Check: isinstance(formation, CompositeFormationContext)? Yes
    ├→ Initialize: total_lbol = 0
    ├→ For each subpop in formation.sub_populations:
    │   ├→ Call: _lbol_atomic(subpop, t)
    │   │   ├→ Call EvolutionContext.lbol(subpop, t)
    │   │   │   ├→ Check subpop type: Simple or Binary?
    │   │   │   ├→ If Simple: use simple_evolution.lbol()
    │   │   │   └→ If Binary: use binary_evolution.lbol()
    │   │   └→ Return subpop luminosity
    │   └→ Add to total_lbol
    └→ Return total_lbol

Note: EvolutionContext automatically handles mixed populations
      (e.g., burst1 = SimpleFormationContext, burst2 = BinaryFormationContext)
```

### Integration Pattern

For population-level properties:

1. **Discrete Mode**: Sum over individual stars
   ```python
   total_lbol = sum(evolution.lbol(m, t) for m in formation.masses)
   ```

2. **Continuous Mode**: Integrate using IMF
   ```python
   total_lbol = integrate(evolution.lbol(m, t) * formation.imf.pdf(m), m_min, m_max)
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
├── population.py                  # StarPopulation, BinaryPopulation (data containers)
├── dist_funcs/                    # IMF and binary distributions
│   ├── imf.py                    # Used by FormationContext
│   └── binaries.py
├── stellar_evolution/            # Used by EvolutionContext
│   ├── isochrone.py
│   ├── se_data_structures.py
│   ├── data_reader.py
│   └── popsynth.py
├── feedbacks/                    # Used by EvolutionContext
│   ├── __init__.py
│   └── sn.py
├── element_yields/               # Used by EvolutionContext
│   ├── yields.py
│   ├── yieldtables.py
│   └── [various implementations]
└── utils/                        # Shared utilities
```

### Dependency Graph

```
SynthPop
├── Formation (one of):
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
└── EvolutionContext
    ├── Evolution Models:
    │   ├── simple_evolution: IsochroneInterpolator (from stellar_evolution/)
    │   └── binary_evolution: BinaryEvolutionModel (from stellar_evolution/) - future
    └── Shared Modules:
        ├── yields: Yields (from element_yields/)
        └── feedback: FBMechanism (from feedbacks/)

Note: CompositeFormationContext does NOT inherit from FormationContext.
      It's a separate container class that holds atomic FormationContext instances.
```

---

## Design Decisions & Rationale

### Why Composition?

- **Modularity**: Contexts can be developed and tested independently
- **Reusability**: Same formation context can be used with different evolution models
- **Extensibility**: New context types (e.g., `BinaryFormationContext`) can be added without modifying core classes

### Why Keep SynthPop as Coordinator?

- **User Convenience**: Simple, high-level API for common use cases
- **Backward Compatibility**: Existing code using `SynthPop` can continue to work
- **Integration Logic**: Population-level calculations naturally belong at this level
- **Context Isolation**: Keeps `FormationContext` and `EvolutionContext` decoupled

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

- Add new `FBMechanism` subclasses for additional feedback types
- Add new `Yields` implementations for different yield tables
- Add new `IMF` classes for different mass functions
- Create custom `EvolutionContext` subclasses for specialized physics

---

## Implementation Status

### Phase 1: Simple Stellar Populations (Current Priority)
- [ ] Create abstract `FormationContext` base class
- [ ] Create `SimpleFormationContext` implementation
- [ ] Create `EvolutionContext` class
- [ ] Refactor `SynthPop` to use composition with `SimpleFormationContext`
- [ ] Implement polymorphic initialization
- [ ] Update tests for new architecture
- [ ] Update documentation and examples
- [ ] Ensure backward compatibility

### Phase 2: Binary Populations (Future)
- [ ] Create `BinaryFormationContext` class
- [ ] Extend `EvolutionContext` for binary evolution
- [ ] Add binary-specific methods to `SynthPop`
- [ ] Integrate binary evolution tracks
- [ ] Add tests for binary populations

### Phase 3: Composite Populations (Future)
- [ ] Create `CompositeFormationContext` class
- [ ] Handle multiple metallicities and ages in `EvolutionContext`
- [ ] Support complex star formation histories
- [ ] Add tests for composite populations

---

## Notes

This architecture document is a living document and will be updated as the implementation progresses. See commit history for evolution of design decisions.
