# Notes on the BPASS Binary Data

These notes record the problems found in the BPASS v2.2.1 stellar models while building
`arsenal_gear.binary_evolution.BPASSReader`, and how the reader handles each one. All numbers are
for the full Z = 0.014 (`z014`) grid: 12,663 binary models and 284 single-star models. Column
numbers are 0-indexed, as in the reader.

The reader uses only the model files (`sneplot-*` in `NEWSINMODS/` and `NEWBINMODS/NEWBINMODS/`).
It does not use BPASS's population input files (`input_bpass_*`, in the IMF tarballs), the
secondary-star models in `NEWBINMODS/NEWSECMODS/`, or the period column. Merger detection,
rejuvenation and the companion's continuation are therefore decided by the reader itself, from the
mass, luminosity, temperature and radius columns.

## Summary

| Issue | Systems affected | Handling in the reader |
|---|---|---|
| Mergers aren't marked | 1,755 mergers (13.9%) | detected from a rise in M1 + M2 (`merger_rise`) |
| Merged primaries gain mass | 402 (3.2%) | kept, flagged in `BinaryTrackSet.mass_excess` |
| Already merged at t = 0 | 245 (1.9%) | merged throughout (`placeholder_mass`) |
| Repeated times | many files | made strictly increasing |
| White-dwarf primaries | 1,775 (14.0%) | companion continues from white-dwarf formation (`wd_log_radius`) |
| Join jumps | 347 jumps > 10% | not handled (#47, #50) |
| "Merger" models without a merger | 224 (1.8%) | used as they are |
| No orbit evolution or eccentricity | all | initial orbit only, e0 = 0 (#45) |

## Mergers aren't marked in the model files

When the primary absorbs the secondary, BPASS adds the secondary's mass to the primary but
leaves the secondary's own columns (mass in col 37, and L, T, R in cols 46–48) unchanged. So
nothing in the file marks the merger, and the secondary appears to keep existing next to the
merged star.

The total mass M1 + M2 is a clean signature, though. A merger shows up as a rise in that total:
- In systems without a merger, the total rises by at most 0.16% in a single step (rounding).
- In mergers, it rises by at least 2.9% in a step, for 99% of them.

The reader treats a single-step rise above `BPASSReader.merger_rise = 0.01` (1%) as a merger. The
rise is relative to the current total. The merger starts at the beginning of the unbroken rise
that contains the first such step; from then on the secondary's mass is 0 and its L, T and R are
NaN. In some mergers (e.g. `5-0.9-0.2`) the secondary's mass column ticks up while it's being
absorbed, so a rule that required the secondary to stay exactly constant would miss them.

BPASS's population input files label models directly (type 0 = merger model, 1 = normal
primary). These labels and the reader's detection agree for all 10,439 type-1 models and for
2,000 of the 2,224 type-0 models; the other 224 are described under "Merger" models without a merger.

## Some merged primaries gain mass

In 402 mergers (3.2% of all systems, 22.9% of mergers; M1 = 1.9–3.7 Msun), BPASS's own merged
primary ends up more massive than the system's initial M1 + M2: by +6.6% to +45%, with a median of
+21% (e.g. `1.9-0.3-0`). This is in BPASS's own primary column, not something the reader adds;
the reader never adds mass.

The values fall into two clean groups. Mergers either stay at or below the initial total (1,353)
or exceed it by at least 6.6% (402), and rounding in the already-merged models (next section)
reaches at most +3.0%. So the 5% threshold `BinaryDataReader.mass_excess_tol` flags exactly these
402. They are kept, and flagged in `BinaryTrackSet.mass_excess`.

## Models already merged at t = 0

245 models (1.9%; M1 = 0.9–2.7 Msun; all type 0) are merged from the first row:
- The primary starts at the combined mass, rounded to BPASS's 0.1 Msun grid, so up to 3.0% above
  the M1 (1 + q) in the model's name.
- The companion is a placeholder: M2 = 1e-5 Msun throughout, with log L = −4 or −25.

The smallest real companion mass in the first row of any other file is 0.01 Msun, so the reader
treats any first-row M2 below `BPASSReader.placeholder_mass = 1e-3` Msun as a placeholder. Such
a system is merged from row 0, with no companion and no single-star continuation. Without this,
the placeholder was continued as a 0.1 Msun star that doesn't exist.

## Repeated and decreasing times

Ages are printed to 7 significant digits, so late phases with short time steps produce repeated
(or slightly decreasing) ages. For example, 944 files have repeated ages during the white-dwarf
cooling phase. The reader makes each model's time strictly increasing with
`utils.array_utils.make_monotonic_increasing`, which spreads the repeated values between their
neighbours. These changes are below the printed precision.

## White-dwarf primaries freeze the companion

When the primary becomes a white dwarf (log R1 < −1), the file goes on through the white dwarf's
cooling, typically for about 100 times longer than the time taken to reach it (a median of 119
times, measured with an earlier merger rule). Over this phase BPASS stops evolving the companion:
- Its columns are completely frozen in 73% of these files.
- In the rest its luminosity drifts only slightly: a median change of 0 dex, and 0.4 dex at the
  99th percentile.
- No companion goes through its own giant phase or death there, even when the cooling phase
  lasts far longer than the companion's lifetime.

While the primary is still alive, BPASS does evolve the companion, approximately: its mass
changes in 70% of files and its luminosity in 84%.

Example: in `3-0.7-3` the primary is a white dwarf at 4.9e8 yr and the file ends at 9.7e10 yr.
Its 2.1 Msun companion, frozen on the main sequence in the file, becomes a giant and then a white
dwarf itself at about 1.3e9 yr on its single-star track.

**Detecting white-dwarf formation:** it's the start of the primary's *final* run of rows below
`BPASSReader.wd_log_radius = -1.0`. The first row below −1 isn't reliable: in 76 systems (mostly
M1 ≈ 2 Msun) the stripped primary dips just below −1, then grows again to log R ≈ 0–1.7 before
finally becoming a white dwarf.

**Handling:**
- The binary phase ends at the white dwarf's formation.
- The primary keeps its whole cooling sequence.
- The companion continues along its single-star track from that point.
- After the join, the rows are the cooling rows plus the single-star rows, with each star
  interpolated onto both.

Of the 1,834 systems without a merger whose primary ends as a white dwarf, 1,775 become one
before the file ends, so their output changes. The other 59 become white dwarfs on the last row.

| Effect on the companion at the join | join at file end (before) | join at white-dwarf formation |
|---|---|---|
| companion's single-star track already over | 381 | 3 |
| companion's mass jump > 50% | 622 | 7 |
| companion's mass jump > 10% | 863 | 347 |
| median jump | 37% | 0.9% |

BPASS's single-star models also continue through white-dwarf cooling, so a companion's
continuation can run to ~1e11 yr. For the same reason, the total ages of low-mass BPASS single
stars aren't comparable with MIST's.

## Mass jumps where the companion joins its single star

When the binary phase ends, the companion continues along the single-star model nearest to its
effective mass:
- normally its initial mass
- if rejuvenated (it accreted more than 5% of its mass and is at least 2 Msun, the BPASS v2.1
  criteria), its maximum mass, as if it was born when it reached that mass

The companion's mass in the binary file and on the single track can differ at the join.

Among the white-dwarf systems:
- 347 jump by more than 10%, and 7 by more than 50%.
- Rejuvenated companions join smoothly: a median jump of 0.6%, and none above 10%.

The larger jumps come from the spacing of the single-star mass grid and from mass transfer after
the companion reached its effective mass. These are modelling choices that could be improved
(#47, #50).

BPASS v2.2 also defines a revised rejuvenation age (its Eq. 1), which the reader doesn't use.

## "Merger" models without a merger

224 models (1.8%) that BPASS labels type 0 ("merger model") show no merger in their files. They
are used as they are, without any special handling.

**What they look like:**
- 222 have M1 ≤ 2 Msun, mostly with q = 0.1–0.3, across the whole period range.
- The companion barely takes part: its mass changes in only 13 files. Where the same M1 and q
  appear at several periods, the primary's track is identical in 14 of 16 cases.
- 212 end at the primary's largest radius:
  - 160 as giants at the tip of the red giant branch (log R1 ≈ 2.2)
  - 54 on or near the main sequence, often at ~1e11 yr
  - 10 as white dwarfs

**Their geometry,** using BPASS's period column (col 34), the current M1 + M2, and Eggleton's
Roche-lobe radius:
- **125 come into contact.** R1 + R2 > a in 100 of them, and the primary overflows its Roche lobe
  in all 125.
  - In 112 of these, the period stays constant (to 0.1%) throughout. The primary grows into a
    ~160 Rsun giant around an orbit as small as 0.02 AU, so (R1 + R2)/a reaches ~40.
  - The other 13 show real interaction through a changing period (e.g. `7-0.9-0.8`: stable mass
    transfer strips the primary to a ~1.4 Msun helium star).
- **99 never come into contact** (M1 = 0.3–1.2 Msun):
  - 47 giants come close (R1/R_L up to 0.98), so later expansion beyond the file's end might
    fill the lobe.
  - 52 end near the main sequence, far inside their lobes (R1/R_L ≤ 0.35).

**Interpretation:** for most contact systems, the file follows the primary as if it were an
isolated star, and BPASS's merger presumably happens outside the detailed model. Some could
instead be contact binaries. This can't be settled from the model files alone; it would need
BPASS's population input files or the papers' description of how they treat these models.

**What the reader does:**
- The primary ends where its file ends (for the giants, at the red-giant tip, with its mass held
  and its L, T and R set to NaN afterwards).
- The companion continues as a single star.
- In the contact systems, the two stars overlap inside an unchanging orbit before that.

## Orbits

The reader stores only each system's initial mass ratio q0, period P0 and eccentricity e0. The
first two are read from the model's name.

BPASS gives the period as it evolves (col 34), but MPA/Bonn files don't, so orbit evolution is
left for later (#45). Neither source gives an eccentricity, so orbits are taken to be circular
(e0 = 0).

## Comparison with MIST single stars

As a check on the single-star models used for the companions' continuations, BPASS and MIST
tracks at 1, 5, 20 and 60 Msun (Z = 0.014 against [Fe/H] = 0) agree well on the main sequence:
- |Δlog L| ≤ 0.009 dex
- |Δlog Teff| ≤ 0.043 dex, with BPASS slightly cooler, most at 60 Msun
- BPASS main sequences last 1.01–1.23 times as long as MIST's

Total lifetimes at 20 and 60 Msun agree to 6%. These comparisons are tested in
`tests/test_be_physics.py`, on the BPASS subset in `tests/data/bpass/`.
