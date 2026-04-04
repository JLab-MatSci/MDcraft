from pathlib import Path

import numpy as np


def read_lammps_data(filename, atoms):
    path = Path(filename)
    if not path.is_file():
        path = Path(__file__).with_name(filename)

    if not path.is_file():
        raise FileNotFoundError(f"LAMMPS data file was not found: {filename}")

    with path.open("r", encoding="utf-8") as stream:
        lines = stream.readlines()

    natoms = None
    bounds = {}
    atoms_start = None

    for idx, raw_line in enumerate(lines):
        line = raw_line.strip()
        if not line:
            continue

        if line.endswith("atoms"):
            natoms = int(line.split()[0])
            continue

        parts = line.split()
        if len(parts) >= 4 and parts[2] in {"xlo", "ylo", "zlo"}:
            lo = float(parts[0])
            hi = float(parts[1])
            axis = parts[2][0]
            bounds[axis] = (lo, hi)
            continue

        if line.startswith("Atoms"):
            atoms_start = idx + 1
            break

    if natoms is None or atoms_start is None or len(bounds) != 3:
        raise ValueError(f"Unsupported or incomplete LAMMPS data file: {path}")

    records = []
    for raw_line in lines[atoms_start:]:
        line = raw_line.strip()
        if not line:
            if records:
                break
            continue
        if line[0].isalpha():
            break

        parts = line.split()
        if len(parts) < 5:
            continue

        records.append((int(parts[0]), float(parts[2]), float(parts[3]), float(parts[4])))
        if len(records) == natoms:
            break

    if len(records) != natoms:
        raise ValueError(f"Expected {natoms} atoms in {path}, found {len(records)}")

    records.sort(key=lambda item: item[0])
    coords = np.asarray([(x, y, z) for _, x, y, z in records], dtype=np.float64)

    atoms.resize(natoms)
    atoms.data["r"]["x"] = coords[:, 0] * 0.1
    atoms.data["r"]["y"] = coords[:, 1] * 0.1
    atoms.data["r"]["z"] = coords[:, 2] * 0.1

    return (
        atoms,
        bounds["x"][0] * 0.1,
        bounds["x"][1] * 0.1,
        bounds["y"][0] * 0.1,
        bounds["y"][1] * 0.1,
        bounds["z"][0] * 0.1,
        bounds["z"][1] * 0.1,
    )
