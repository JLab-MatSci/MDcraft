"""MPI-safe thermodynamic reductions and root-only global LAMMPS dump."""

from pathlib import Path

import numpy as np
from mpi4py import MPI

KB = 0.008314462175  # kJ / (mol*K); original analysis.py
ATOMIC_MASS_GRAM = 1.6605390666e-24
NA = 6.02214076e23


def global_thermodynamics(atoms, volume, comm):
    """Return global (Ep,Et,Ek,T,virial,P,density) using only owned atoms.

    Formulas / units match the MPI analysis example provided by the user.
    Ghost/alien atoms must never be included in reductions.
    """
    data = atoms.data
    mass = np.asarray(data["m"], dtype=np.float64)
    vx = np.asarray(data["v"]["x"], dtype=np.float64)
    vy = np.asarray(data["v"]["y"], dtype=np.float64)
    vz = np.asarray(data["v"]["z"], dtype=np.float64)
    momentum_local = np.array([np.sum(mass * vx), np.sum(mass * vy),
                               np.sum(mass * vz), np.sum(mass)], dtype=np.float64)
    momentum_global = np.zeros(4, dtype=np.float64)
    comm.Allreduce(momentum_local, momentum_global, op=MPI.SUM)
    if momentum_global[3] <= 0:
        raise ValueError("global atom mass is zero")
    center = momentum_global[:3] / momentum_global[3]
    two_ek = np.sum(mass * ((vx - center[0]) ** 2 +
                            (vy - center[1]) ** 2 + (vz - center[2]) ** 2))
    energy = float(np.sum(data["Ep"]))
    if "Em" in data.dtype.names:
        energy += float(np.sum(data["Em"]))
    trace = float(np.sum(data["V"]["xx"] + data["V"]["yy"] + data["V"]["zz"]))
    local = np.array([two_ek, energy, trace, np.sum(mass), mass.size], dtype=np.float64)
    global_stats = np.zeros(5, dtype=np.float64)
    comm.Allreduce(local, global_stats, op=MPI.SUM)
    two_ek, ep, virial, total_mass, n_atoms = global_stats
    volume_cm3 = volume * 1e-21
    if n_atoms <= 0 or volume_cm3 <= 0:
        raise ValueError("empty atom system or invalid domain volume")
    ek = two_ek / 2.0
    temperature = two_ek / (3.0 * KB * n_atoms)
    pressure = (two_ek - virial) / (3.0 * volume_cm3 * NA)
    density = ATOMIC_MASS_GRAM * total_mass / volume_cm3
    return ep, ep + ek, ek, temperature, virial, pressure, density


SNAPSHOT_DTYPE = np.dtype([
    ("id", "u4"), ("type", "i4"),
    ("x", "f8"), ("y", "f8"), ("z", "f8"),
])


def local_snapshot(atoms):
    """Create compact owned-atom array from structured MDCraft atoms."""
    data = atoms.data
    snapshot = np.empty(data.size, dtype=SNAPSHOT_DTYPE)
    snapshot["id"] = data["uid"]
    snapshot["type"] = np.asarray(data["kind"], dtype=np.int32) + 1  # Al=1, Ni=2
    for axis in ("x", "y", "z"):
        snapshot[axis] = data["r"][axis]
    if not np.all((snapshot["type"] == 1) | (snapshot["type"] == 2)):
        raise ValueError("unexpected atomic kind (expected Al=0, Ni=1)")
    return snapshot


def combine_snapshots(chunks, expected_atoms=None):
    data = np.concatenate(chunks) if chunks else np.empty(0, dtype=SNAPSHOT_DTYPE)
    data.sort(order="id")
    if expected_atoms is not None and len(data) != expected_atoms:
        raise ValueError(f"gather produced {len(data)} atoms; expected {expected_atoms}")
    if data.size and (data["id"][0] < 1 or np.any(data["id"][1:] == data["id"][:-1])):
        raise ValueError("duplicate or invalid atom UIDs in MPI snapshot")
    return data


def write_lammps_snapshot(path, data, bounds, step, coordinate_scale=10.0):
    """Write ordinary LAMMPS text dump. Angstrom scale is explicit, not implicit."""
    xmin, xmax, ymin, ymax, zmin, zmax = bounds
    scale = coordinate_scale
    with Path(path).open("w", encoding="utf-8") as output:
        output.write(f"ITEM: TIMESTEP\n{step}\nITEM: NUMBER OF ATOMS\n{len(data)}\n")
        output.write("ITEM: BOX BOUNDS pp pp pp\n")
        for low, high in ((xmin, xmax), (ymin, ymax), (zmin, zmax)):
            output.write(f"{low * scale:.14g} {high * scale:.14g}\n")
        output.write("ITEM: ATOMS id type x y z\n")
        if len(data):
            values = np.column_stack((data["id"], data["type"], data["x"] * scale,
                                      data["y"] * scale, data["z"] * scale))
            np.savetxt(output, values, fmt=["%.0f", "%.0f", "%.12g", "%.12g", "%.12g"])


def dump_frame(atoms, comm, bounds, dump_dir, step, expected_atoms, coordinate_scale=10.0):
    """Collect owned atoms from all ranks; only root writes a complete snapshot."""
    # Check per-rank errors BEFORE gather, and root I/O errors AFTER gather:
    # otherwise other MPI ranks could hang indefinitely at the next collective.
    error = None
    try:
        local = local_snapshot(atoms)
    except Exception as exc:
        error = f"rank {comm.Get_rank()}: {exc}"
        local = None
    errors = comm.allgather(error)
    if any(errors):
        raise RuntimeError("invalid local snapshot: " + "; ".join(e for e in errors if e))
    chunks = comm.gather(local, root=0)
    error = None
    if comm.Get_rank() == 0:
        try:
            full = combine_snapshots(chunks, expected_atoms)
            path = Path(dump_dir) / f"frame{step:09d}.lammpstrj"
            write_lammps_snapshot(path, full, bounds, step, coordinate_scale)
        except Exception as exc:
            error = f"rank 0 output error: {exc}"
    error = comm.bcast(error, root=0)
    if error:
        raise RuntimeError(error)
