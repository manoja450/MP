"""Helpers for IAEA-format phase-space files (e.g. Varian TrueBeam phase spaces).

IAEA phase spaces come as a pair ``NAME.IAEAheader`` + ``NAME.IAEAphsp``.
TOPAS expects ``NAME.header`` + ``NAME.phsp`` and is pointed at ``NAME``
(without extension), so :func:`link_for_topas` creates matching symlinks.

The header is a plain-text file organised in ``$SECTION:`` blocks, e.g.::

    $RECORD_CONTENTS:
        1     // X is stored ?
        1     // Y is stored ?
        0     // Z is stored ?
    $RECORD_CONSTANT:
        26.7000   // Constant Z
    $ORIG_HISTORIES:
        123456789

When Z is *not* stored, every particle sits on the same plane and TOPAS places
it at z = 0 of the source component, so we must position that component at the
plane given by the ``Constant Z`` value (distance from the target, in cm).
"""
from __future__ import annotations

import os
import re
from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class PhspInfo:
    path_stem: str
    record_length: int | None = None
    orig_histories: int | None = None
    particles: int | None = None
    photons: int | None = None
    electrons: int | None = None
    positrons: int | None = None
    z_stored: bool = True
    weight_stored: bool = True
    z_constant_cm: float | None = None
    sections: dict[str, list[str]] = field(default_factory=dict, repr=False)

    def summary(self) -> str:
        lines = [f"Phase space: {self.path_stem}"]
        for k in ("record_length", "orig_histories", "particles", "photons", "electrons", "positrons"):
            lines.append(f"  {k:16s}: {getattr(self, k)}")
        lines.append(f"  Z stored        : {self.z_stored}")
        lines.append(f"  Weight stored   : {self.weight_stored}")
        lines.append(f"  Constant Z (cm) : {self.z_constant_cm}")
        return "\n".join(lines)


_SECTION_RE = re.compile(r"^\s*\$([A-Z0-9_]+):\s*(.*)$")


def _split_sections(text: str) -> dict[str, list[str]]:
    sections: dict[str, list[str]] = {}
    current = None
    for raw in text.splitlines():
        m = _SECTION_RE.match(raw)
        if m:
            current = m.group(1)
            sections[current] = []
            if m.group(2).strip():
                sections[current].append(m.group(2).strip())
            continue
        if current is not None and raw.strip():
            sections[current].append(raw.strip())
    return sections


def _first_number(lines: list[str]) -> float | None:
    for line in lines:
        m = re.match(r"^\s*([-+]?\d+(?:\.\d*)?(?:[eE][-+]?\d+)?)", line)
        if m:
            return float(m.group(1))
    return None


def _value_and_comment(line: str) -> tuple[str, str]:
    if "//" in line:
        v, c = line.split("//", 1)
        return v.strip(), c.strip()
    return line.strip(), ""


def find_header(stem: str | os.PathLike) -> Path:
    stem = str(stem)
    for ext in (".IAEAheader", ".header"):
        p = Path(stem + ext)
        if p.exists():
            return p
    raise FileNotFoundError(f"No .IAEAheader or .header found for {stem}")


def read_header(stem: str | os.PathLike) -> PhspInfo:
    """Parse the IAEA header belonging to ``stem`` (path without extension)."""
    hdr = find_header(stem)
    sec = _split_sections(hdr.read_text(errors="replace"))
    info = PhspInfo(path_stem=str(stem), sections=sec)

    def num(name):
        v = _first_number(sec.get(name, []))
        return None if v is None else int(v)

    info.record_length = num("RECORD_LENGTH")
    info.orig_histories = num("ORIG_HISTORIES")
    info.particles = num("PARTICLES")
    info.photons = num("PHOTONS")
    info.electrons = num("ELECTRONS")
    info.positrons = num("POSITRONS")

    for line in sec.get("RECORD_CONTENTS", []):
        val, com = _value_and_comment(line)
        low = com.lower()
        if low.startswith("z is stored"):
            info.z_stored = val.startswith("1")
        elif low.startswith("weight is stored"):
            info.weight_stored = val.startswith("1")

    for line in sec.get("RECORD_CONSTANT", []):
        val, com = _value_and_comment(line)
        if re.search(r"\bconstant\s+z\b", com.lower()):
            info.z_constant_cm = float(val)
    return info


def link_for_topas(stem: str | os.PathLike, dest_dir: str | os.PathLike) -> Path:
    """Create ``dest_dir/NAME.header`` and ``dest_dir/NAME.phsp`` symlinks.

    Returns the new stem (``dest_dir/NAME``) to give to TOPAS.
    """
    stem = Path(stem)
    dest_dir = Path(dest_dir)
    dest_dir.mkdir(parents=True, exist_ok=True)
    pairs = [(".IAEAheader", ".header"), (".IAEAphsp", ".phsp")]
    for src_ext, dst_ext in pairs:
        src = Path(str(stem) + src_ext)
        if not src.exists():
            src = Path(str(stem) + dst_ext)  # already TOPAS-style
        if not src.exists():
            raise FileNotFoundError(f"Missing {stem}{src_ext}")
        dst = dest_dir / (stem.name + dst_ext)
        if dst.is_symlink() or dst.exists():
            dst.unlink()
        dst.symlink_to(src.resolve())
    return dest_dir / stem.name


# ---------------------------------------------------------------------------
# Splitting a phase space into independent chunks
# ---------------------------------------------------------------------------
# Record layout (IAEA format, as read by OpenTOPAS TsSourcePhaseSpace):
#   byte 0      : int8 particle type (1 photon, 2 electron, 3 positron, ...);
#                 its sign encodes the sign of the z direction cosine
#   bytes 1..4  : float32 energy; a negative value flags a NEW HISTORY
#   then x, y, [z], u, v, [weight], extra floats, extra longs
# Only the first 5 bytes are inspected here; records are copied verbatim.

import json

import numpy as np

_COUNT_SECTIONS = {"PARTICLES": "particles", "PHOTONS": "photons",
                   "ELECTRONS": "electrons", "POSITRONS": "positrons",
                   "ORIG_HISTORIES": "orig_histories"}


def _data_file(stem: str | os.PathLike) -> Path:
    for ext in (".IAEAphsp", ".phsp"):
        p = Path(str(stem) + ext)
        if p.exists():
            return p
    raise FileNotFoundError(f"No .IAEAphsp or .phsp found for {stem}")


def count_records(records: np.ndarray) -> dict[str, int]:
    """Count particle types and new-history flags in an (n, record_length) uint8 array."""
    ptype = np.abs(records[:, 0].view(np.int8).astype(np.int16))
    energy = records[:, 1:5].copy().view("<f4").ravel()
    return {
        "particles": int(records.shape[0]),
        "photons": int(np.count_nonzero(ptype == 1)),
        "electrons": int(np.count_nonzero(ptype == 2)),
        "positrons": int(np.count_nonzero(ptype == 3)),
        "new_history_flags": int(np.count_nonzero(energy < 0)),
    }


def _rewrite_header(text: str, values: dict[str, int]) -> str:
    """Replace the first number after each $SECTION: line with a new value."""
    out, pending = [], None
    for line in text.splitlines():
        m = _SECTION_RE.match(line)
        if m:
            pending = m.group(1) if m.group(1) in values else None
            out.append(line)
            continue
        if pending and line.strip():
            num = re.match(r"^(\s*)([-+]?\d+(?:\.\d*)?(?:[eE][-+]?\d+)?)(.*)$", line)
            if num:
                line = f"{num.group(1)}{values[pending]}{num.group(3)}"
                pending = None
        out.append(line)
    return "\n".join(out) + "\n"


def split_phase_space(stem: str | os.PathLike, out_dir: str | os.PathLike, n_chunks: int,
                      block_records: int = 2_000_000) -> list[dict]:
    """Split one IAEA phase space into ``n_chunks`` contiguous, TOPAS-readable chunks.

    Writes ``out_dir/<name>_cNNN.header`` and ``.phsp`` plus a ``chunks.json``
    manifest. ``histories`` in the manifest is the number of new-history flags in
    the chunk or, if the file has none (Varian files), the number of photons, which
    matches TOPAS' ``LimitedAssumePhotonIsNewHistory`` behaviour.
    """
    info = read_header(stem)
    if info.record_length is None:
        raise ValueError("header has no $RECORD_LENGTH")
    order = " ".join(info.sections.get("BYTE_ORDER", [])).strip()
    if order and not order.startswith("1234"):
        raise ValueError(f"only little-endian (1234) files are supported, header says {order!r}")
    rl = info.record_length
    data = _data_file(stem)
    size = data.stat().st_size
    if size % rl:
        raise ValueError(f"file size {size} is not a multiple of record length {rl}")
    n_total = size // rl
    if n_chunks < 1 or n_chunks > n_total:
        raise ValueError("n_chunks must be between 1 and the number of records")
    header_text = find_header(stem).read_text(errors="replace")
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    mm = np.memmap(data, dtype=np.uint8, mode="r", shape=(n_total, rl))
    edges = np.linspace(0, n_total, n_chunks + 1).astype(np.int64)
    name = Path(str(stem)).name
    manifest = []
    for c in range(n_chunks):
        a, b = int(edges[c]), int(edges[c + 1])
        counts = {"particles": 0, "photons": 0, "electrons": 0, "positrons": 0, "new_history_flags": 0}
        chunk_stem = out_dir / f"{name}_c{c:03d}"
        with open(str(chunk_stem) + ".phsp", "wb") as fout:
            for s in range(a, b, block_records):
                block = np.asarray(mm[s:min(b, s + block_records)])
                for k, v in count_records(block).items():
                    counts[k] += v
                fout.write(block.tobytes())
        hist = counts["new_history_flags"] or counts["photons"]
        orig = info.orig_histories
        orig_chunk = int(round(orig * (b - a) / n_total)) if orig else hist
        new_vals = {"PARTICLES": counts["particles"], "PHOTONS": counts["photons"],
                    "ELECTRONS": counts["electrons"], "POSITRONS": counts["positrons"],
                    "ORIG_HISTORIES": orig_chunk, "CHECKSUM": (b - a) * rl}
        Path(str(chunk_stem) + ".header").write_text(_rewrite_header(header_text, new_vals))
        manifest.append({"stem": str(chunk_stem), "source": str(stem), "first_record": a,
                         "n_records": b - a, "histories": int(hist),
                         "orig_histories_scaled": orig_chunk, **counts})
    with open(out_dir / f"{name}_chunks.json", "w") as f:
        json.dump({"source": str(stem), "record_length": rl, "z_stored": info.z_stored,
                   "z_constant_cm": info.z_constant_cm, "chunks": manifest}, f, indent=2)
    return manifest


def load_chunk_manifests(chunk_dir: str | os.PathLike) -> list[dict]:
    """Read every ``*_chunks.json`` in ``chunk_dir`` and return all chunks (sorted)."""
    chunks = []
    meta = None
    for p in sorted(Path(chunk_dir).glob("*_chunks.json")):
        d = json.loads(p.read_text())
        meta = meta or d
        if (d["z_stored"], d["z_constant_cm"]) != (meta["z_stored"], meta["z_constant_cm"]):
            raise ValueError(f"{p}: phase-space plane differs from other files")
        for ch in d["chunks"]:
            ch["z_stored"] = d["z_stored"]
            ch["z_constant_cm"] = d["z_constant_cm"]
            chunks.append(ch)
    if not chunks:
        raise FileNotFoundError(f"no *_chunks.json manifests in {chunk_dir}; run scripts/split_phsp.py first")
    return chunks


# ---------------------------------------------------------------------------
# Sampling particles (to check coordinates before trusting a phase space)
# ---------------------------------------------------------------------------
_STORED_ORDER = ["X", "Y", "Z", "U", "V", "W", "Weight"]


def record_layout(info: PhspInfo) -> list[str]:
    """Names of the float32 fields after the type byte and energy, in file order."""
    stored = {k: True for k in _STORED_ORDER}
    stored["W"] = False                       # W is normally derived from U, V
    n_extra_f = n_extra_l = 0
    for line in info.sections.get("RECORD_CONTENTS", []):
        val, com = _value_and_comment(line)
        low = com.lower()
        for k in _STORED_ORDER:
            if k != "W" and low.startswith(f"{k.lower()} is stored"):
                stored[k] = val.startswith("1")   # W is never in the record (as OpenTOPAS reads it)
        if low.startswith("extra floats stored"):
            n_extra_f = int(float(val))
        if low.startswith("extra longs stored"):
            n_extra_l = int(float(val))
    fields = [k for k in _STORED_ORDER if stored[k]]
    fields += [f"extra_float_{i}" for i in range(n_extra_f)]
    fields += [f"extra_long_{i}" for i in range(n_extra_l)]
    return fields


def sample_particles(stem, n: int = 1_000_000) -> dict[str, np.ndarray]:
    """Read the first ``n`` particles into a dict of arrays (cm, MeV)."""
    info = read_header(stem)
    fields = record_layout(info)
    expected = 1 + 4 + 4 * len(fields)
    if expected != info.record_length:
        raise ValueError(f"layout {fields} gives {expected} bytes but header says {info.record_length}")
    data = _data_file(stem)
    n = min(n, data.stat().st_size // info.record_length)
    raw = np.fromfile(data, dtype=np.uint8, count=n * info.record_length).reshape(n, info.record_length)
    ptype = raw[:, 0].view(np.int8).astype(int)
    floats = raw[:, 1:].copy().view("<f4").reshape(n, -1)
    out = {"type": np.abs(ptype), "cosz_negative": ptype < 0,
           "energy": np.abs(floats[:, 0]), "new_history": floats[:, 0] < 0}
    for i, f in enumerate(fields):
        out[f] = floats[:, 1 + i]
    return out


def describe_sample(s: dict) -> str:
    def rng(a):
        return f"min {a.min():9.3f}  mean {a.mean():9.3f}  max {a.max():9.3f}"
    lines = [f"particles sampled      : {s['energy'].size}"]
    for k, name in ((1, "photons"), (2, "electrons"), (3, "positrons")):
        lines.append(f"  {name:9s}            : {np.mean(s['type'] == k) * 100:6.2f} %")
    lines.append(f"energy (MeV)           : {rng(s['energy'])}")
    for k in ("X", "Y", "Z", "U", "V"):
        if k in s:
            lines.append(f"{k:22s} : {rng(s[k])}")
    if "Z" in s:
        z = s["Z"]
        lines.append(f"distinct Z values      : {np.unique(np.round(z, 4)).size} (first: {np.unique(np.round(z, 4))[:5]})")
    lines.append(f"sign(cos z) negative   : {np.mean(s['cosz_negative']) * 100:6.2f} % of particles")
    lines.append(f"new-history flags      : {np.mean(s['new_history']) * 100:6.2f} %")
    if "U" in s and "V" in s:
        r = np.hypot(s["X"], s["Y"]) if "X" in s else None
        w = np.sqrt(np.clip(1 - s["U"] ** 2 - s["V"] ** 2, 0, 1))
        if r is not None:
            sel = r > 5
            radial = (s["X"][sel] * s["U"][sel] + s["Y"][sel] * s["V"][sel]) / r[sel]
            lines.append(f"mean radial direction  : {radial.mean():+.4f} (positive = beam diverging outward)")
        lines.append(f"mean |cos z|           : {w.mean():.4f}")
    return "\n".join(lines)
