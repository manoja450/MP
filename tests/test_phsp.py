import json

import numpy as np
import pytest

from mcverify.phsp import (count_records, link_for_topas, load_chunk_manifests, read_header,
                           split_phase_space)
from tests.helpers import write_fake_phsp


def test_read_header(tmp_path):
    info = write_fake_phsp(tmp_path / "ps" / "TB6X_00", n=500)
    h = read_header(tmp_path / "ps" / "TB6X_00")
    assert h.record_length == info["rl"]
    assert h.particles == 500
    assert h.photons == info["photons"]
    assert h.z_stored is False
    assert h.weight_stored is True
    assert h.z_constant_cm == pytest.approx(26.7)


def test_split_preserves_every_record(tmp_path):
    stem = tmp_path / "ps" / "TB6X_00"
    info = write_fake_phsp(stem, n=1003)
    chunks = split_phase_space(stem, tmp_path / "chunks", n_chunks=4, block_records=100)
    assert len(chunks) == 4
    assert sum(c["n_records"] for c in chunks) == 1003
    assert sum(c["photons"] for c in chunks) == info["photons"]
    joined = b"".join((tmp_path / "chunks" / f"TB6X_00_c{c:03d}.phsp").read_bytes() for c in range(4))
    assert joined == (tmp_path / "ps" / "TB6X_00.IAEAphsp").read_bytes()
    # rewritten headers are consistent with their chunk
    h = read_header(tmp_path / "chunks" / "TB6X_00_c001")
    assert h.particles == chunks[1]["n_records"]
    assert h.photons == chunks[1]["photons"]
    assert h.z_constant_cm == pytest.approx(26.7)
    # no new-history flags in the file -> histories = photons
    assert chunks[0]["histories"] == chunks[0]["photons"]


def test_new_history_flags_are_counted(tmp_path):
    stem = tmp_path / "TB"
    info = write_fake_phsp(stem, n=300, new_history_flags=True)
    chunks = split_phase_space(stem, tmp_path / "c", n_chunks=1)
    assert chunks[0]["new_history_flags"] == info["flags"]
    assert chunks[0]["histories"] == info["flags"]


def test_manifest_loading(tmp_path):
    for k in range(2):
        write_fake_phsp(tmp_path / f"TB_{k:02d}", n=200, seed=k)
        split_phase_space(tmp_path / f"TB_{k:02d}", tmp_path / "chunks", n_chunks=3)
    chunks = load_chunk_manifests(tmp_path / "chunks")
    assert len(chunks) == 6
    assert all(c["z_constant_cm"] == pytest.approx(26.7) for c in chunks)
    assert json.loads((tmp_path / "chunks" / "TB_00_chunks.json").read_text())["record_length"] == 25


def test_link_for_topas(tmp_path):
    write_fake_phsp(tmp_path / "src" / "TB", n=10)
    new = link_for_topas(tmp_path / "src" / "TB", tmp_path / "linked")
    assert (tmp_path / "linked" / "TB.header").is_symlink()
    assert (tmp_path / "linked" / "TB.phsp").read_bytes() == (tmp_path / "src" / "TB.IAEAphsp").read_bytes()
    assert new == tmp_path / "linked" / "TB"


def test_count_records_types():
    rec = np.zeros((4, 25), dtype=np.uint8)
    rec[:, 0] = np.array([1, -1, 2, -3], dtype=np.int8).view(np.uint8)
    e = np.array([1.0, -2.0, 3.0, 4.0], dtype="<f4")
    rec[:, 1:5] = e.view(np.uint8).reshape(4, 4)
    c = count_records(rec)
    assert (c["photons"], c["electrons"], c["positrons"], c["new_history_flags"]) == (2, 1, 1, 1)


def test_sample_particles_layout(tmp_path):
    from mcverify.phsp import describe_sample, record_layout, sample_particles
    write_fake_phsp(tmp_path / "TB", n=500)
    info = read_header(tmp_path / "TB")
    assert record_layout(info) == ["X", "Y", "U", "V", "Weight"]
    s = sample_particles(tmp_path / "TB", n=100)
    assert s["energy"].size == 100 and "Z" not in s
    assert (s["energy"] > 0).all()
    assert "photons" in describe_sample(s)


def test_sample_particles_with_z(tmp_path):
    from mcverify.phsp import record_layout, sample_particles
    stem = tmp_path / "TBZ"
    n, rl = 50, 25
    rec = np.zeros((n, rl), dtype=np.uint8)
    rec[:, 0] = np.full(n, 1, dtype=np.int8).view(np.uint8)
    fl = np.column_stack([np.full(n, 2.0), np.zeros(n), np.zeros(n), np.full(n, 26.7),
                          np.zeros(n), np.zeros(n)]).astype("<f4")
    rec[:, 1:] = fl.view(np.uint8).reshape(n, 24)
    (tmp_path / "TBZ.IAEAphsp").write_bytes(rec.tobytes())
    hdr = ("$RECORD_CONTENTS:\n 1 // X is stored ?\n 1 // Y is stored ?\n 1 // Z is stored ?\n"
           " 1 // U is stored ?\n 1 // V is stored ?\n 1 // W is stored ?\n 0 // Weight is stored ?\n"
           " 0 // Extra floats stored ?\n 0 // Extra longs stored ?\n$RECORD_LENGTH:\n 25\n$PARTICLES:\n 50\n")
    (tmp_path / "TBZ.IAEAheader").write_text(hdr)
    assert record_layout(read_header(stem)) == ["X", "Y", "Z", "U", "V"]
    s = sample_particles(stem)
    assert np.allclose(s["Z"], 26.7) and np.allclose(s["energy"], 2.0)


def test_manifest_with_relative_paths_is_made_absolute(tmp_path, monkeypatch):
    write_fake_phsp(tmp_path / "TB", n=100)
    monkeypatch.chdir(tmp_path)
    split_phase_space("TB", "rel_chunks", n_chunks=2)       # relative paths, as typed by a user
    monkeypatch.chdir("/")
    chunks = load_chunk_manifests(tmp_path / "rel_chunks")
    assert all(c["stem"].startswith(str(tmp_path)) for c in chunks)
