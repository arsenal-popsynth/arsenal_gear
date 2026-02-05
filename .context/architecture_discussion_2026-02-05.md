# Architecture Discussion - February 5, 2026

## Summary

This conversation focused on designing the architecture for the Arsenal Gear codebase refactoring. We created a comprehensive architecture document at [docs/architecture_claude.md](../docs/architecture_claude.md) that specifies the proposed design.

## Key Decisions Made

### 1. Top-Level Architecture (Composition Pattern)

**Decision**: Use composition over inheritance with three main classes:
- **`SynthPop`** - Main user-facing API (renamed from `StellarPopulation`)
- **`FormationContext`** - Encapsulates initial conditions and birth properties
- **`EvolutionContext`** - Encapsulates time-dependent evolution and physical processes

**Rationale**:
- Clear separation of concerns
- Flexibility to swap contexts independently
- Easier testing of individual components
- Ability to reuse contexts across multiple populations

### 2. SynthPop Public Interface

Defined comprehensive public interface organized by physical process:

**Supernova Methods**:
- `nsn(t)` - Cumulative supernovae
- `ndotsn(t)` - Supernova rate

**Luminosity Methods**:
- Population integrated: `lbol(t)`, `lum(t, λ_min, λ_max)`, `mag(t, filter)`
- Per-star (discrete mode): `lbol_stars(t)`, `lum_stars(t, λ_min, λ_max)`, `mag_stars(t, filter)`

**Mechanical Feedback Methods**:
- Population integrated: `edot_wind(t)`, `edot_sn(t)`, `edot_total(t)`
- Per-star (discrete mode): `edot_wind_stars(t)`, `mdot_wind_stars(t)`

**Mass Loss Methods**:
- `mdot_wind(t)`, `mloss_wind(t)`, `mloss_sn(t)`, `mloss_total(t)`

**Chemical Yields Methods**:
- `yields(element, t)`, `yields_by_source(element, t)`

**Stellar Property Methods**:
- Population integrated: `teff(t)`
- Per-star (discrete mode): `teff_stars(t)`, `mass_stars(t)`

**Note**: Removed `mmax(t)` from public interface - it's an internal implementation detail.

### 3. FormationContext Architecture (Container Pattern)

**Decision**: Use abstract base class for atomic populations + separate container for composite populations.

**Architecture Pattern**:
```
FormationContext (abstract base - for atomic populations only)
├── SimpleFormationContext (single IMF, metallicity, age)
└── BinaryFormationContext (future - includes binary parameters)

CompositeFormationContext (container - does NOT inherit from FormationContext)
└── Contains List[FormationContext] (atomic contexts only)
```

**Key Design Choice**: Atomic contexts return scalar values for metallicity/age, while composite returns arrays. This avoids the problem of atomic contexts having to return arrays of identical values.

**Rationale**:
1. **Conceptual Clarity**: An atomic population *has* a metallicity (not an array)
2. **Type Correctness**: Scalars for atomic, arrays for composite
3. **Natural Composition**: Composite literally contains other contexts
4. **Easy Querying**: "What's the metallicity?" has clear meaning
5. **Straightforward Iteration**: `SynthPop` loops over sub-populations for composite

**How SynthPop Handles Both**:
- Type checks: `isinstance(formation, CompositeFormationContext)`
- For atomic: computes directly
- For composite: iterates over `sub_populations` and sums results

### 4. EvolutionContext Architecture (Type Dispatch)

**Decision**: Single `EvolutionContext` class with internal type dispatch.

**Structure**:
```python
class EvolutionContext:
    # Type-specific evolution models
    simple_evolution: IsochroneInterpolator      # For SimpleFormationContext
    binary_evolution: BinaryEvolutionModel       # For BinaryFormationContext (future)

    # Shared modules (work for both)
    yields: Yields
    feedback: FBMechanism
```

**Dispatch Pattern**:
```python
def lbol(self, formation: FormationContext, t):
    if isinstance(formation, SimpleFormationContext):
        return self._lbol_simple(formation, t)
    elif isinstance(formation, BinaryFormationContext):
        return self._lbol_binary(formation, t)
```

**Alternative Considered**: Parallel `SimpleEvolutionContext` and `BinaryEvolutionContext` subclasses.

**Rationale for Single Class**:
1. **User-friendly**: One context "just works" for all population types
2. **Works with CompositeFormationContext**: Handles mixed simple/binary populations
3. **Shared modules**: Yields and feedback reused for both types
4. **Consistent pattern**: Matches how `SynthPop` dispatches on atomic vs composite
5. **Existing code structure**: Modules already separated (IsochroneInterpolator, yields, feedback)

### 5. Naming Convention

**Decision**: Rename `StellarPopulation` → `SynthPop`

**Rationale**: Makes clearer distinction between:
- `SynthPop` = Top-level population synthesis coordinator
- Stellar populations = Individual atomic formation contexts

## Implementation Strategy

### Phase 1: Simple Stellar Populations (Current Priority)
- [ ] Create abstract `FormationContext` base class
- [ ] Create `SimpleFormationContext` implementation
- [ ] Create `EvolutionContext` class (with `simple_evolution` module)
- [ ] Refactor `SynthPop` to use composition
- [ ] Implement polymorphic initialization
- [ ] Update tests
- [ ] Ensure backward compatibility

### Phase 2: Binary Populations (Future)
- [ ] Create `BinaryFormationContext` class
- [ ] Add `binary_evolution` module to `EvolutionContext`
- [ ] Add binary-specific methods to `SynthPop`
- [ ] Integrate binary evolution tracks

### Phase 3: Composite Populations (Future)
- [ ] Create `CompositeFormationContext` container class
- [ ] Support complex star formation histories
- [ ] Add tests for composite populations

## File Organization

```
arsenal_gear/
├── __init__.py                    # SynthPop (main API)
├── formation_context/             # NEW
│   ├── __init__.py               # Abstract FormationContext base
│   ├── simple.py                 # SimpleFormationContext
│   ├── binary.py                 # BinaryFormationContext (future)
│   └── composite.py              # CompositeFormationContext (future)
├── evolution_context.py           # EvolutionContext class (NEW)
├── population.py                  # StarPopulation, BinaryPopulation (data containers)
├── dist_funcs/                    # IMF and binary distributions
├── stellar_evolution/             # Used by EvolutionContext
├── feedbacks/                     # Used by EvolutionContext
├── element_yields/                # Used by EvolutionContext
└── utils/                         # Shared utilities
```

## Key Architectural Patterns

### Pattern 1: Composition Over Inheritance
- `SynthPop` composes `FormationContext` and `EvolutionContext`
- Not using inheritance hierarchy

### Pattern 2: Container Pattern (FormationContext)
- Atomic contexts (`SimpleFormationContext`, `BinaryFormationContext`) inherit from abstract base
- `CompositeFormationContext` is separate container, does NOT inherit
- Container holds list of atomic contexts

### Pattern 3: Type Dispatch (EvolutionContext)
- Single class with multiple evolution models
- Dispatches internally based on formation context type
- Shared modules (yields, feedback) for all types

### Pattern 4: Polymorphic Initialization (SynthPop)
```python
# Option 1: Simple parameters (auto-constructs contexts)
pop = SynthPop(Mtot=1e6, metallicity=0.0, ...)

# Option 2: Explicit contexts (power users)
formation = SimpleFormationContext(...)
evolution = EvolutionContext(...)
pop = SynthPop(formation_context=formation, evolution_context=evolution)

# Option 3: Mixed
formation = SimpleFormationContext(...)
pop = SynthPop(formation_context=formation, metallicity=0.0, ...)
```

## Discussion Points and Rationale

### Why Container Pattern for CompositeFormationContext?

**Question**: Should `CompositeFormationContext` return arrays for metallicity/age, or should it be a container?

**Decision**: Container pattern - it holds multiple atomic contexts with their own scalar metallicity/age values.

**Pros of Container Approach**:
- Conceptually clearer: each sub-population has definite metallicity/age
- Easy to query individual sub-populations
- Simpler atomic contexts (no array wrapping needed)
- Type checking is straightforward

**Cons**:
- Requires `isinstance()` checks in `SynthPop`
- Different code paths for atomic vs composite

**Conclusion**: Conceptual clarity and ease of use outweigh the minor complexity of type checking.

### Why Single EvolutionContext for Simple and Binary?

**Question**: Should we have separate `SimpleEvolutionContext` and `BinaryEvolutionContext` classes?

**Decision**: Single `EvolutionContext` with internal dispatch.

**Pros of Single Class**:
- User creates one context that works for everything
- Naturally handles `CompositeFormationContext` with mixed populations
- Shared yields/feedback modules - no duplication
- Matches the dispatch pattern used elsewhere

**Cons**:
- Loads both evolution models even if only using one
- Type checking required internally

**Conclusion**: User-friendliness and consistency with composite pattern outweigh minor overhead.

## Next Steps

1. Start implementing Phase 1 (Simple Stellar Populations)
2. Begin with `FormationContext` abstract base class
3. Implement `SimpleFormationContext`
4. Create `EvolutionContext` with `simple_evolution` module
5. Refactor existing `StellarPopulation` → `SynthPop` using new architecture

## References

- Main architecture document: [docs/architecture_claude.md](../docs/architecture_claude.md)
- Codebase summary: [.context/CODEBASE_SUMMARY.md](CODEBASE_SUMMARY.md)
- Current implementation: [arsenal_gear/__init__.py](../arsenal_gear/__init__.py)

## Context for Future Conversations

When resuming this work:
1. Review [docs/architecture_claude.md](../docs/architecture_claude.md) for the full specification
2. Check [.context/CODEBASE_SUMMARY.md](CODEBASE_SUMMARY.md) for current code structure
3. The architecture is fully specified but **not yet implemented**
4. Phase 1 (Simple Stellar Populations) is the current priority
5. All design decisions have been documented with rationale
