# Arsenal Gear Codebase Summary

**Last Updated**: 2026-02-05

## Quick Reference

**Project**: Arsenal Gear - Lightweight population synthesis code for stellar populations with emphasis on feedback from massive stars

**Main Entry Point**: [arsenal_gear/__init__.py](arsenal_gear/__init__.py) - Contains `StellarPopulation` class

**Status**: Partial refactoring in progress. Core functionality working, but needs architectural refactoring to implement composition pattern described in [docs/architecture_claude.md](docs/architecture_claude.md)

---

## Current Architecture vs. Proposed

### Current (Monolithic)
```
StellarPopulation
├─ Handles IMF (hardcoded Salpeter)
├─ Samples masses
├─ Stores Mtot, metallicity
├─ Manages IsochroneInterpolator
└─ Computes population properties
```

### Proposed (Composition)
```
StellarPopulation (coordinator)
├─ FormationContext
│  ├─ IMF (pluggable)
│  ├─ Sampled masses
│  ├─ Mtot, metallicity
│  └─ Binary parameters
└─ EvolutionContext
   ├─ IsochroneInterpolator
   ├─ Yields
   └─ FBMechanism
```

**Key Difference**: Proposed architecture separates formation (initial conditions) from evolution (time-dependent physics) using composition pattern.

---

## Directory Structure

```
arsenal_gear/
├── arsenal_gear/
│   ├── __init__.py                    # StellarPopulation (main API)
│   ├── population.py                  # StarPopulation, BinaryPopulation (data containers)
│   ├── dist_funcs/                    # IMF and binary distributions
│   │   ├── imf.py                     # IMF base, Salpeter implementation
│   │   └── binaries.py                # Binary distributions
│   ├── stellar_evolution/             # Isochrone interpolation
│   │   ├── isochrone.py               # IsochroneInterpolator
│   │   ├── se_data_structures.py      # Isochrone, Track data structures
│   │   └── data_reader.py             # MIST data reader
│   ├── feedbacks/                     # Feedback mechanisms
│   │   ├── __init__.py                # FBMechanism base class
│   │   └── sn.py                      # Supernova functions
│   ├── element_yields/                # Chemical yields
│   │   ├── yields.py                  # Yields abstract base
│   │   ├── yieldtables.py             # YieldTables base
│   │   ├── source.py                  # Yield interpolation
│   │   ├── limongichieffi2018.py      # LC18 yield tables
│   │   └── nugrid.py                  # NuGrid yields
│   └── utils/                         # Utilities
├── tests/                             # Test suite
├── docs/                              # Sphinx documentation
│   └── architecture_claude.md         # Proposed architecture doc
└── input/                             # Example parameter files
```

---

## Key Classes and Locations

### Core Classes

| Class | File | Status | Purpose |
|-------|------|--------|---------|
| `StellarPopulation` | [arsenal_gear/__init__.py](arsenal_gear/__init__.py) | ✅ Working | Main user-facing API |
| `StarPopulation` | [arsenal_gear/population.py](arsenal_gear/population.py#L10) | ✅ Working | Star data container |
| `BinaryPopulation` | [arsenal_gear/population.py](arsenal_gear/population.py#L50) | ⚠️ Not integrated | Binary data container |

### Formation/Sampling

| Class | File | Status | Purpose |
|-------|------|--------|---------|
| `IMF` | [arsenal_gear/dist_funcs/imf.py](arsenal_gear/dist_funcs/imf.py) | ✅ Working | Base IMF class |
| `Salpeter` | [arsenal_gear/dist_funcs/imf.py](arsenal_gear/dist_funcs/imf.py) | ✅ Working | Salpeter IMF implementation |
| `BinaryDistribution` | [arsenal_gear/dist_funcs/binaries.py](arsenal_gear/dist_funcs/binaries.py) | ⚠️ Not integrated | Binary property sampling |

### Stellar Evolution

| Class | File | Status | Purpose |
|-------|------|--------|---------|
| `IsochroneInterpolator` | [arsenal_gear/stellar_evolution/isochrone.py](arsenal_gear/stellar_evolution/isochrone.py) | ✅ Working | Interpolate stellar properties |
| `Isochrone` | [arsenal_gear/stellar_evolution/se_data_structures.py](arsenal_gear/stellar_evolution/se_data_structures.py) | ✅ Working | Isochrone data structure |
| `IsochroneSet` | [arsenal_gear/stellar_evolution/se_data_structures.py](arsenal_gear/stellar_evolution/se_data_structures.py) | ✅ Working | Collection of isochrones |
| `StellarTrack` | [arsenal_gear/stellar_evolution/se_data_structures.py](arsenal_gear/stellar_evolution/se_data_structures.py) | ✅ Working | Single stellar track |
| `TrackSet` | [arsenal_gear/stellar_evolution/se_data_structures.py](arsenal_gear/stellar_evolution/se_data_structures.py) | ✅ Working | Collection of tracks |
| `MISTReader` | [arsenal_gear/stellar_evolution/data_reader.py](arsenal_gear/stellar_evolution/data_reader.py) | ✅ Working | Load MIST data |

### Yields

| Class | File | Status | Purpose |
|-------|------|--------|---------|
| `Yields` | [arsenal_gear/element_yields/yields.py](arsenal_gear/element_yields/yields.py) | ⚠️ Not integrated | Abstract yield interface |
| `YieldTables` | [arsenal_gear/element_yields/yieldtables.py](arsenal_gear/element_yields/yieldtables.py) | ⚠️ Not integrated | Base for tabulated yields |
| `LimongiChieffi2018` | [arsenal_gear/element_yields/limongichieffi2018.py](arsenal_gear/element_yields/limongichieffi2018.py) | ⚠️ Not integrated | LC18 yield tables |
| `Source` | [arsenal_gear/element_yields/source.py](arsenal_gear/element_yields/source.py) | ⚠️ Not integrated | Multi-dimensional yield interpolation |

### Feedback

| Class | File | Status | Purpose |
|-------|------|--------|---------|
| `FBMechanism` | [arsenal_gear/feedbacks/__init__.py](arsenal_gear/feedbacks/__init__.py) | ❌ Interface only | Abstract feedback base |
| `lifetimes_Raiteri` | [arsenal_gear/feedbacks/sn.py](arsenal_gear/feedbacks/sn.py) | ✅ Working | Stellar lifetimes function |

---

## Current StellarPopulation Class

**Location**: [arsenal_gear/__init__.py](arsenal_gear/__init__.py)

### Key Attributes
```python
Mtot              # Total population mass (Quantity)
metallicity       # log10(Z/Zsun)
discrete          # Boolean: discrete vs continuous mode
imf               # IMF object (hardcoded Salpeter)
masses            # Sampled stellar masses (if discrete)
Nstar             # Number of stars
iso               # IsochroneInterpolator instance
verbose           # Debug flag
```

### Public Methods
- `lbol(t)` - Bolometric luminosity at time t
- `teff(t)` - Luminosity-weighted effective temperature
- `nsn(t)` - Cumulative number of supernovae
- `ndotsn(t)` - Supernova rate
- `lbol_iso(t)` - Per-star bolometric luminosity
- `teff_iso(t)` - Per-star effective temperature
- `_integrate_pop(iso, q)` - Internal: integrate quantity over population

### Current Limitations
1. IMF is **hardcoded** to Salpeter (0.08-100 Msun, α=2.3)
2. Metallicity is **fixed at initialization**
3. No **yield integration** (yields classes exist but not used)
4. No **feedback integration** (FBMechanism only interface)
5. **Binary populations** not supported (BinaryPopulation class unused)
6. Tight coupling between formation and evolution logic

---

## Important Patterns and Conventions

### Units
- Heavy use of **Astropy Quantity** objects throughout
- Always use `.to(u.unit)` for conversions
- Dimensionless: `u.dimensionless_unscaled`

### Naming Conventions
- `log_*` or `l*` prefix for logarithmic quantities (e.g., `lteff_name`, `llbol`)
- `mini` = initial mass
- `Mmax` = maximum mass of living stars
- `Lbol` = bolometric luminosity
- `Teff` = effective temperature
- `eep` = Equivalent Evolutionary Phase

### Masking
- Uses `numpy.ma.masked_array` for out-of-bounds values
- Masks applied when masses outside interpolation range

### Data Containers
- `StarPopulation` and `BinaryPopulation` inherit from `dict`
- Use dataclasses for structured data (`Isochrone`, `StellarTrack`)

---

## Data Flow Examples

### Initialization
```
User creates StellarPopulation(Mtot=1e6, metallicity=0.0, discrete=True)
    ↓
__init__() creates:
    ├─ Salpeter IMF (hardcoded: 0.08-100 Msun, α=2.3)
    ├─ Samples masses if discrete=True
    └─ IsochroneInterpolator(met=metallicity)
        └─ MISTReader loads MIST data
```

### Query: `pop.lbol(t)`
```
pop.lbol(t)
    ├─ Discrete mode:
    │  └─ Sum: iso.lbol(mass, t) for each mass in self.masses
    │
    └─ Continuous mode:
       └─ Integrate: ∫ iso.lbol(m, t) * imf.pdf(m) dm
```

### Query: `pop.nsn(t)`
```
pop.nsn(t)
    ├─ Get Mmax(t) from iso.mmax(t)
    ├─ Compute fraction of stars > 8 Msun (hardcoded threshold)
    └─ Return Nstar * fraction
```

---

## What Needs to Change for Proposed Architecture

### 1. Create `FormationContext` Class
**New file**: [arsenal_gear/formation_context.py](arsenal_gear/formation_context.py)

**Responsibilities**:
- Store Mtot, metallicity, IMF
- Sample stellar masses (discrete mode)
- Provide mass distributions
- Handle binary properties (future)

**Attributes**:
```python
Mtot: Quantity["mass"]
metallicity: float
imf: IMF
masses: Quantity["mass"]  # if discrete
Nstar: int
discrete: bool
```

**Methods**:
```python
sample_masses() -> Quantity["mass"]
get_mass_distribution() -> tuple
expected_stars_above(mass) -> float
```

### 2. Create `EvolutionContext` Class
**New file**: [arsenal_gear/evolution_context.py](arsenal_gear/evolution_context.py)

**Responsibilities**:
- Coordinate IsochroneInterpolator, Yields, FBMechanism
- Provide time-dependent stellar properties
- Compute yields and feedback

**Composed Components**:
```python
stellar_evolution: IsochroneInterpolator
yields: Yields
feedback: FBMechanism
```

**Methods**:
```python
lbol(masses, metallicity, t) -> Quantity["power"]
teff(masses, metallicity, t) -> Quantity["temperature"]
mmax(t) -> Quantity["mass"]
sn_rate(masses, t) -> float
yields_at_time(element, masses, metallicity, t) -> Quantity["mass"]
```

### 3. Refactor `StellarPopulation`
**Changes to**: [arsenal_gear/__init__.py](arsenal_gear/__init__.py)

**New initialization pattern**:
```python
# Option 1: Simple (auto-construct contexts)
pop = StellarPopulation(Mtot=1e6, metallicity=0.0, ...)

# Option 2: Explicit contexts
formation = FormationContext(...)
evolution = EvolutionContext(...)
pop = StellarPopulation(formation_context=formation, evolution_context=evolution)

# Option 3: Mixed
pop = StellarPopulation(formation_context=formation, metallicity=0.0, ...)
```

**New attributes**:
```python
formation: FormationContext
evolution: EvolutionContext
```

**Updated methods** - delegate to contexts:
```python
def lbol(self, t):
    masses = self.formation.masses  # or mass distribution
    return integrate(self.evolution.lbol(masses, self.formation.metallicity, t))
```

### 4. Integrate Yields and Feedback
- Connect `Yields` classes to `EvolutionContext`
- Implement concrete `FBMechanism` subclasses
- Add methods to `StellarPopulation`:
  - `yields(element, t)`
  - `feedback_energy(t)`

### 5. Update Tests
- Add tests for `FormationContext`
- Add tests for `EvolutionContext`
- Update `StellarPopulation` tests for new initialization
- Add integration tests

---

## Testing

**Test files**:
- [tests/test_constructors.py](tests/test_constructors.py) - StarPopulation, BinaryPopulation
- [tests/test_imf.py](tests/test_imf.py) - IMF sampling
- [tests/test_binaries.py](tests/test_binaries.py) - Binary distributions
- [tests/test_isochrone.py](tests/test_isochrone.py) - Isochrone interpolation

**Test patterns**:
- Use `numpy.testing.assert_array_equal` for array comparisons
- `IsochroneInterpolator(test=True)` for cross-validation
- Compare discrete vs continuous mode results

---

## Implementation Roadmap

From [docs/architecture_claude.md](docs/architecture_claude.md):

- [ ] Create `FormationContext` class
- [ ] Create `EvolutionContext` class
- [ ] Refactor `StellarPopulation` to use composition
- [ ] Implement polymorphic initialization
- [ ] Integrate yields into `EvolutionContext`
- [ ] Integrate feedback into `EvolutionContext`
- [ ] Update tests for new architecture
- [ ] Update documentation and examples
- [ ] Ensure backward compatibility

---

## Key Implementation Details

### Isochrone Interpolation
- Two modes: `interp_op="iso"` (isochrone-based) or `"track"` (EEP-based)
- Uses **EEP (Equivalent Evolutionary Phase)** to match evolutionary stages
- PCHIP interpolation preserves monotonicity
- Supports MIST stellar evolution models

### IMF Sampling
- Inherits from `scipy.stats.rv_continuous`
- `sample_mass(Mtot)` estimates N ≈ Mtot / mean_mass
- Returns Astropy Quantity with Msun units

### Supernova Calculation
- Uses `iso.mmax(t)` to find maximum living mass
- Hardcoded **8 Msun** minimum for CCSN (should be configurable)
- `nsn(t)` integrates IMF above Mmax
- `ndotsn(t)` uses chain rule with dMmax/dt

---

## Common Gotchas

1. **IMF is hardcoded** - Cannot currently configure IMF in StellarPopulation
2. **8 Msun threshold** - CCSN minimum mass is hardcoded, not configurable
3. **Metallicity is fixed** - Cannot vary with time
4. **Yields not integrated** - Yield classes exist but aren't used by StellarPopulation
5. **No feedback physics** - FBMechanism is only an interface
6. **Binary populations unused** - BinaryPopulation class exists but not integrated

---

## Quick Start Guide

### Current Usage
```python
from arsenal_gear import StellarPopulation
import astropy.units as u

# Create a population
pop = StellarPopulation(
    Mtot=1e6 * u.Msun,
    metallicity=0.0,  # Solar metallicity
    discrete=True,    # Sample discrete stars
    verbose=False
)

# Query properties
import numpy as np
times = np.logspace(6, 10, 100) * u.yr

lbol = pop.lbol(times)      # Bolometric luminosity
teff = pop.teff(times)      # Effective temperature
nsn = pop.nsn(times)        # Cumulative supernovae
sn_rate = pop.ndotsn(times) # Supernova rate
```

### Proposed Usage (After Refactor)
```python
from arsenal_gear import StellarPopulation, FormationContext, EvolutionContext
from arsenal_gear.dist_funcs.imf import Kroupa
from arsenal_gear.element_yields import LimongiChieffi2018
import astropy.units as u

# Option 1: Simple
pop = StellarPopulation(Mtot=1e6*u.Msun, metallicity=0.0)

# Option 2: Custom IMF
formation = FormationContext(Mtot=1e6*u.Msun, imf=Kroupa(...))
pop = StellarPopulation(formation_context=formation)

# Option 3: Custom yields
yields = LimongiChieffi2018()
evolution = EvolutionContext(yields=yields)
pop = StellarPopulation(formation_context=formation, evolution_context=evolution)

# New methods
pop.yields('O', times)      # Oxygen yields
pop.yields('Fe', times)     # Iron yields
pop.feedback_energy(times)  # Feedback energy
```

---

## Notes

- This summary is based on codebase exploration on 2026-02-05
- Architecture proposal documented in [docs/architecture_claude.md](docs/architecture_claude.md)
- Current implementation is functional but needs refactoring for proposed architecture
- Git status shows untracked files: `arsenal_gear/element_yields/data/`, `docs/architecture_claude.md`, `playground.ipynb`
- Main branch: `main`
