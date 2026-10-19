"""
make_bpass_fixture.py
=====================

Regenerates the small BPASS subset in tests/data/bpass/ used by
tests/test_be_physics.py, from a local copy of the BPASS v2.2.1 models
(the extracted NEWSINMODS/ and NEWBINMODS/ folders).

Usage:
    python tests/data/make_bpass_fixture.py [BPASS_ROOTDIR]

BPASS_ROOTDIR defaults to arsenal_gear/binary_evolution/data/BPASS.
"""

import shutil
import sys
from pathlib import Path

import numpy as np

from arsenal_gear.binary_evolution import BPASSReader

MET = "z014"
# single stars compared with MIST
SINGLE_MASSES = ["1", "5", "20", "60"]
# range of primary masses the example binaries are drawn from
M1_RANGE = (10, 30)
# also included: a merger whose secondary column ticks up while it is absorbed,
# a merger whose primary ends up more massive than the initial M1 + M2, and a
# model already merged at t = 0 (placeholder companion)
EXTRA_BINARIES = [
    "sneplot-z014-5-0.9-0.2",
    "sneplot-z014-1.9-0.3-0",
    "sneplot-z014-1.1-0.5-0",
]
# columns kept: the reader uses columns 1-5 of single stars and 1-48 of binaries
SINGLE_COLS, BINARY_COLS = 6, 49

HERE = Path(__file__).parent
OUT = HERE / "bpass"


def trim(src: Path, dst: Path, ncols: int) -> None:
    """Copy the first ncols columns of a model file, as text (values unchanged)."""
    dst.parent.mkdir(parents=True, exist_ok=True)
    with open(src, encoding="utf-8") as fin, open(dst, "w", encoding="utf-8") as fout:
        for line in fin:
            fout.write(" ".join(line.split()[:ncols]) + "\n")


def nearest(masses: dict, mass: float) -> str:
    """Name of the single star whose mass is nearest to mass."""
    return min(masses, key=lambda name: abs(masses[name] - mass))


def classify(path: Path) -> tuple[str, float]:
    """
    Kind of binary ("merger", "rejuvenated" or "plain"), using the reader's
    criteria, and the companion's effective mass (for its continuation).
    """
    p_mass, s_mass = np.genfromtxt(path, usecols=(5, 37)).T
    total = p_mass + s_mass
    if np.any(np.diff(total) / total[:-1] > BPASSReader.merger_rise):
        return "merger", np.nan
    m2_init, m2_max = s_mass[0], s_mass.max()
    if (m2_max > 1.05 * m2_init) and (m2_max >= 2):
        return "rejuvenated", m2_max
    return "plain", m2_init


def main(rootdir: Path) -> None:
    single_dir = rootdir / "NEWSINMODS" / MET
    binary_dir = rootdir / "NEWBINMODS" / "NEWBINMODS" / MET
    single_files = {f.name.split("-")[-1]: f for f in single_dir.glob("sneplot-*")}
    single_masses = {name: float(name) for name in single_files}

    chosen = {}
    for path in sorted(binary_dir.glob("sneplot-*")):
        m1 = float(path.name.split("-")[-3])
        if not M1_RANGE[0] <= m1 <= M1_RANGE[1]:
            continue
        kind, m2_eff = classify(path)
        if kind not in chosen:
            chosen[kind] = (path, m2_eff)
        if len(chosen) == 3:
            break

    singles = set(SINGLE_MASSES)
    for kind, (path, m2_eff) in chosen.items():
        if kind != "merger":
            # the single star the companion continues onto, as the reader picks it
            singles.add(nearest(single_masses, m2_eff))
        print(f"{kind}: {path.name}")

    # replace the model folders, keeping the README
    for folder in ("NEWSINMODS", "NEWBINMODS"):
        shutil.rmtree(OUT / folder, ignore_errors=True)
    for name in sorted(singles, key=float):
        trim(single_files[name], OUT / "NEWSINMODS" / MET / single_files[name].name,
             SINGLE_COLS)
    binaries = [path for path, _ in chosen.values()]
    binaries += [binary_dir / name for name in EXTRA_BINARIES]
    for path in binaries:
        trim(path, OUT / "NEWBINMODS" / "NEWBINMODS" / MET / path.name, BINARY_COLS)
    size = sum(f.stat().st_size for f in OUT.rglob("*") if f.is_file())
    print(f"single stars: {sorted(singles, key=float)}; total size {size / 1e6:.2f} MB")


if __name__ == "__main__":
    default = HERE.parents[1] / "arsenal_gear" / "binary_evolution" / "data" / "BPASS"
    main(Path(sys.argv[1]) if len(sys.argv) > 1 else default)
