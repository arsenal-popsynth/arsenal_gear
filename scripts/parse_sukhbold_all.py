#!/usr/bin/env python3
"""
Sukhbold 2016 Full Yield + Energy Parser (corrected version)
===========================================================

Parses:
- Nucleosynthesis yields from *<mass>.yield_table files
  in N20/Z9.6 subdirs
- Explosion energies from the results_N20 / results_Z9.6 files

Outputs three CSVs:
- sukhbold_ccsn_yields.csv      (ejecta = CCSN contribution)
- sukhbold_wind_yields.csv      (pre-SN wind contribution)
- sukhbold_explosion_energies.csv

Default paths:
  Yields:  ../data/Sukhbold2016/nucleosynthesis_yields/{N20,Z9.6}/[sz]<mass>.yield_table
  Energies: ../data/Sukhbold2016/explosion_results_PHOTB/results_{N20,Z9.6}

Usage:
    python parse_sukhbold_all.py
"""

import argparse
import re
from pathlib import Path

import pandas as pd


def parse_yield_table(filepath: Path):
    """
    Parse one *<mass>.yield_table file
    Returns mass, dict of isotope → (ejecta, wind)
    """
    text = filepath.read_text(encoding="utf-8").strip()

    # Extract mass from filename (s12.25.yield_table or z15.yield_table → 12.25 or 15)
    mass_match = re.search(r"[sz](\d+\.?\d*)\.yield_table", filepath.name)
    if not mass_match:
        print(f"Skipping {filepath} — no mass in filename")
        return None, None, None

    mass = float(mass_match.group(1))

    ejecta = {}
    wind = {}

    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("["):
            continue
        parts = re.split(r"\s+", line)
        if len(parts) >= 3:
            isotope = parts[0].lower()
            try:
                ej = float(parts[1])
                wi = float(parts[2])
                ejecta[isotope] = ej
                wind[isotope] = wi
            except ValueError:
                continue

    return mass, ejecta, wind


def parse_energy_file(filepath: Path):
    """
    Parse the big results_N20 / results_Z9.6 file for energies.
    """
    text = filepath.read_text(encoding="utf-8")

    data = []
    model_matches = list(
        re.finditer(
            r"^([a-zA-Z]?\d+\.?\d*)\s*(?:\([^\)]*\))?\s*-*$", text, re.MULTILINE
        )
    )

    for i, match in enumerate(model_matches):
        start = match.start()
        end = model_matches[i + 1].start() if i + 1 < len(model_matches) else len(text)
        block = text[start:end]

        mass_match = re.search(r"\d+\.?\d*", match.group(1))
        if not mass_match:
            continue
        mass = float(mass_match.group())

        e_match = re.search(r"E_exp\s*=\s*([\d\.]+)\s*foe", block, re.IGNORECASE)
        if e_match:
            energy = float(e_match.group(1))
            data.append({"mass_Msun": mass, "explosion_energy": energy})

    if not data:
        print(f"No energies parsed from {filepath}")
        return pd.DataFrame()

    return pd.DataFrame(data).sort_values("mass_Msun").reset_index(drop=True)


def collect_yields_from_dir(yield_dir: Path):
    """
    Find all s*.yield_table or z*.yield_table files in a directory and parse them.
    Returns two wide DataFrames: ejecta and wind
    """
    files = sorted(yield_dir.glob("[sz]*.yield_table"))
    if not files:
        print(f"No [sz]*.yield_table files found in {yield_dir}")
        return pd.DataFrame(), pd.DataFrame()

    ejecta_list = []
    wind_list = []

    for f in files:
        mass, ej_dict, wi_dict = parse_yield_table(f)
        if mass is None:
            continue

        ej_row = {"mass_Msun": mass}
        wi_row = {"mass_Msun": mass}

        ej_row.update({f"{k}_ejecta": v for k, v in ej_dict.items()})
        wi_row.update({f"{k}_wind": v for k, v in wi_dict.items()})

        ejecta_list.append(ej_row)
        wind_list.append(wi_row)

    df_e = pd.DataFrame(ejecta_list).sort_values("mass_Msun").reset_index(drop=True)
    df_w = pd.DataFrame(wind_list).sort_values("mass_Msun").reset_index(drop=True)

    return df_e, df_w


def main():
    default_base = Path(__file__).resolve().parent.parent / "data" / "Sukhbold2016"

    parser = argparse.ArgumentParser(
        description="Parse all Sukhbold 2016 yields and energies"
    )
    parser.add_argument(
        "--n20-yields-dir",
        default=default_base / "nucleosynthesis_yields" / "N20",
        type=Path,
    )
    parser.add_argument(
        "--z96-yields-dir",
        default=default_base / "nucleosynthesis_yields" / "Z9.6",
        type=Path,
    )
    parser.add_argument(
        "--n20-energy-file",
        default=default_base / "explosion_results_PHOTB" / "results_N20",
        type=Path,
    )
    parser.add_argument(
        "--z96-energy-file",
        default=default_base / "explosion_results_PHOTB" / "results_Z9.6",
        type=Path,
    )
    parser.add_argument("--output-dir", default=default_base, type=Path)

    args = parser.parse_args()

    print(f"N20 yields dir:   {args.n20_yields_dir}")
    print(f"Z9.6 yields dir:  {args.z96_yields_dir}")
    print(f"N20 energy file:  {args.n20_energy_file}")
    print(f"Z9.6 energy file: {args.z96_energy_file}")
    print(f"Output dir:       {args.output_dir}")

    # Collect yields
    print("\nCollecting nucleosynthesis yields...")
    df_e_n20, df_w_n20 = collect_yields_from_dir(args.n20_yields_dir)
    df_e_z96, df_w_z96 = collect_yields_from_dir(args.z96_yields_dir)

    # Parse energies
    print("Parsing explosion energies...")
    en_n20 = parse_energy_file(args.n20_energy_file)
    en_z96 = parse_energy_file(args.z96_energy_file)

    # Combine
    df_ccsn = (
        pd.concat([df_e_z96, df_e_n20]).sort_values("mass_Msun").reset_index(drop=True)
    )
    df_wind = (
        pd.concat([df_w_z96, df_w_n20]).sort_values("mass_Msun").reset_index(drop=True)
    )
    df_energy = (
        pd.concat([en_z96, en_n20]).sort_values("mass_Msun").reset_index(drop=True)
    )

    # Save
    args.output_dir.mkdir(parents=True, exist_ok=True)
    df_ccsn.to_csv(args.output_dir / "sukhbold_ccsn_yields.csv", index=False)
    df_wind.to_csv(args.output_dir / "sukhbold_wind_yields.csv", index=False)
    df_energy.to_csv(args.output_dir / "sukhbold_explosion_energies.csv", index=False)

    print("\nDone. Saved:")
    print(f"  CCSN yields:   {len(df_ccsn)} rows")
    print(f"  Wind yields:   {len(df_wind)} rows")
    print(f"  Energies:      {len(df_energy)} rows")


if __name__ == "__main__":
    main()
