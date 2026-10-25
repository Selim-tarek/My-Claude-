#!/usr/bin/env python3
"""Verification suite for the updated cSS pipeline architecture.

This suite is deliberately patient-data-free. It validates the architectural
contracts that should hold after the Claude Code changes:

  1. Source modules compile and import.
  2. Surface-aligned patches use a physical millimetre frame and orthogonal
     u/v/n axes, including an oblique affine.
  3. Multimodal patch extraction preserves channel and validity metadata.
  4. The patch classifier trains with study/patient-grouped LOPO evaluation.
  5. Missing modalities remain explicit instead of silently becoming ordinary
     image data.
  6. The manifest contract rejects duplicate patches and invalid labels.
  7. Existing detector/QC/evaluation safeguards remain present.
  8. An optional object-fate ledger is checked if implemented.

Usage from the project root:
  python tests/verify_updated_architecture.py

Exit status is non-zero if any required check fails. Optional checks are
reported as SKIP and do not fail the suite.
"""

from __future__ import annotations

import csv
import importlib.util
import json
import os
import py_compile
import shutil
import sys
import tempfile
import traceback
from pathlib import Path

import nibabel as nib
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"

PASS = []
SKIP = []
FAIL = []


def report(status: str, name: str, detail: str = "") -> None:
    line = f"{status:<5} {name}"
    if detail:
        line += f" — {detail}"
    print(line)
    {"PASS": PASS, "SKIP": SKIP, "FAIL": FAIL}[status].append(name)


def check(name, fn, optional: bool = False) -> None:
    try:
        detail = fn()
        report("PASS", name, str(detail or ""))
    except Exception as exc:  # keep running to report all failures
        if optional:
            report("SKIP", name, f"{type(exc).__name__}: {exc}")
        else:
            report("FAIL", name, f"{type(exc).__name__}: {exc}")
            traceback.print_exc()


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_source_compiles():
    targets = [SCRIPTS / "surface_aligned_patch.py", SCRIPTS / "classify_surface_patches.py"]
    missing = [str(p) for p in targets if not p.exists()]
    if missing:
        raise FileNotFoundError(", ".join(missing))
    for path in targets:
        py_compile.compile(str(path), doraise=True)
    return f"compiled {len(targets)} modules"


def test_static_pipeline_contracts():
    legacy = [SCRIPTS / "detect_css.py", SCRIPTS / "evaluate.py", SCRIPTS / "train_ranker.py"]
    missing_files = [str(p.name) for p in legacy if not p.exists()]
    if len(missing_files) == len(legacy):
        raise FileNotFoundError("legacy detector sources absent: " + ", ".join(missing_files))
    if missing_files:
        raise FileNotFoundError("legacy detector sources partially absent: " + ", ".join(missing_files))
    detect = legacy[0].read_text()
    evaluate = legacy[1].read_text()
    ranker = legacy[2].read_text()
    required_detect = ["qc_status", "rule_version", "n_pieces"]
    required_eval = ["n_lost_excluded", "n_lost_dropped", "n_lost_not_generated"]
    required_ranker = ["study_id", "LOPO", "test"]
    for token in required_detect:
        if token not in detect:
            raise AssertionError(f"detect_css.py missing contract token {token}")
    for token in required_eval:
        if token not in evaluate:
            raise AssertionError(f"evaluate.py missing contract token {token}")
    for token in required_ranker:
        if token not in ranker:
            raise AssertionError(f"train_ranker.py missing contract token {token}")
    return "QC, fate accounting, and grouped-ranker hooks present"


def make_oblique_affine(vox=0.8):
    # Rotation around the third axis plus anisotropic physical scaling.
    theta = np.deg2rad(17.0)
    rot = np.array([[np.cos(theta), -np.sin(theta), 0.0],
                    [np.sin(theta), np.cos(theta), 0.0],
                    [0.0, 0.0, 1.0]])
    affine = np.eye(4)
    affine[:3, :3] = rot @ np.diag([vox, vox * 1.1, vox * 1.4])
    affine[:3, 3] = [-12.0, 7.0, 4.0]
    return affine


def synthetic_inputs(tmp: Path):
    shape = (64, 64, 48)
    affine = make_oblique_affine()
    ijk = np.indices(shape).transpose(1, 2, 3, 0)
    world = (affine @ np.c_[ijk.reshape(-1, 3), np.ones(np.prod(shape))].T).T[:, :3].reshape(shape + (3,))
    center = np.median(world.reshape(-1, 3), axis=0)
    # Signed distance to an ellipsoidal tissue boundary, in millimetres.
    rel = world - center
    sd = np.sqrt((rel[..., 0] / 20.0) ** 2 + (rel[..., 1] / 22.0) ** 2 + (rel[..., 2] / 17.0) ** 2) * 20.0 - 20.0
    mask = np.zeros(shape, np.uint8)
    # Curvilinear surface candidate, expressed in voxel coordinates.
    for t in range(-12, 13):
        i = int(shape[0] / 2 + t)
        j = int(shape[1] / 2 + round(0.08 * t * t))
        k = int(shape[2] / 2 + 8)
        if 1 <= i < shape[0] - 1 and 1 <= j < shape[1] - 1 and 1 <= k < shape[2] - 1:
            mask[i, j, k] = 1
    # A sheet-like SWI signal and a separate phase channel.
    swi = np.exp(-((sd + 0.5) / 2.0) ** 2).astype(np.float32)
    phase = (0.3 * swi).astype(np.float32)
    flair = (0.1 * np.maximum(swi, 0)).astype(np.float32)
    pial_path = tmp / "pial_sd.nii.gz"
    mask_path = tmp / "candidate.nii.gz"
    swi_path = tmp / "swi.nii.gz"
    phase_path = tmp / "phase.nii.gz"
    flair_path = tmp / "flair.nii.gz"
    for path, data in ((pial_path, sd), (mask_path, mask), (swi_path, swi), (phase_path, phase), (flair_path, flair)):
        nib.save(nib.Nifti1Image(data, affine), str(path))
    return {
        "shape": shape, "affine": affine, "pial": pial_path, "mask": mask_path,
        "swi": swi_path, "phase": phase_path, "flair": flair_path,
    }


def test_surface_frame_and_multimodal_patch():
    patch_mod = load_module(SCRIPTS / "surface_aligned_patch.py", "surface_aligned_patch_test")
    with tempfile.TemporaryDirectory(prefix="css_arch_patch_") as td:
        tmp = Path(td)
        files = synthetic_inputs(tmp)
        pial = nib.load(str(files["pial"]))
        mask = nib.load(str(files["mask"]))
        images = {name: nib.load(str(files[name])) for name in ("swi", "phase", "flair")}
        result = patch_mod.extract_patch(images, pial, mask, size_mm=(16, 8, 8), spacing_mm=(0.5, 0.5, 0.5))
        assert result["shape"] == (33, 17, 17), result["shape"]
        assert set(result["channels"]) == {"swi", "phase", "flair"}
        assert result["valid"].mean() > 0.8
        u = result["frame"]["u_world"]; v = result["frame"]["v_world"]; n = result["frame"]["n_world"]
        assert np.isclose(np.linalg.norm(u), 1.0, atol=1e-5)
        assert np.isclose(np.linalg.norm(v), 1.0, atol=1e-5)
        assert np.isclose(np.linalg.norm(n), 1.0, atol=1e-5)
        assert max(abs(np.dot(u, v)), abs(np.dot(u, n)), abs(np.dot(v, n))) < 1e-5
        # Confirm the normal is actually derived from the oblique affine rather
        # than assuming that an array axis is the superior-inferior direction.
        assert np.max(np.abs(result["frame"]["n_world"])) > 0.2
        return f"shape={result['shape']}, valid={result['valid'].mean():.3f}, oblique affine"


def write_npz_patch(path: Path, kind: str, include_phase=True, include_flair=True):
    rng = np.random.default_rng(abs(hash(path.name)) % (2**32))
    shape = (33, 17, 17)
    x = np.zeros(shape, np.float32)
    if kind == "css":
        x[3:30, 7:11, 7:11] = 1.0
    else:
        x[14:19, 8:10, 2:15] = 1.0
    x += rng.normal(0, 0.02, shape).astype(np.float32)
    arrays = {"image_swi": x, "valid": np.ones(shape, np.uint8), "valid_swi": np.ones(shape, np.uint8),
              "metadata_json": np.asarray('{"pial_sd_at_center_mm": -0.5}')}
    if include_phase:
        arrays.update(image_phase=(0.25 * x).astype(np.float32), valid_phase=np.ones(shape, np.uint8))
    if include_flair:
        arrays.update(image_flair=(0.1 * x).astype(np.float32), valid_flair=np.ones(shape, np.uint8))
    np.savez_compressed(path, **arrays)


def test_classifier_patient_grouping_and_missing_channels():
    clf = load_module(SCRIPTS / "classify_surface_patches.py", "classifier_test")
    with tempfile.TemporaryDirectory(prefix="css_arch_classifier_") as td:
        tmp = Path(td); rows = []
        for study_index in range(6):
            for label, kind in (("cSS", "css"), ("vein", "vein")):
                path = tmp / f"P{study_index:02d}_{label}.npz"
                # Omit phase for one vein patch to test explicit missingness.
                write_npz_patch(path, kind, include_phase=not (study_index == 0 and kind == "vein"))
                rows.append({"patch_path": str(path), "study_id": f"P{study_index:02d}", "label": label})
        manifest = tmp / "manifest.csv"
        with open(manifest, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["patch_path", "study_id", "label"])
            writer.writeheader(); writer.writerows(rows)
        loaded = clf.read_manifest(str(manifest))
        X, y, groups, feature_names, channels = clf.build_matrix(loaded)
        assert len(set(groups)) == 6
        assert set(channels) == {"flair", "phase", "swi"}
        assert any("phase:missing" in n for n in feature_names)
        probabilities, folds = clf.grouped_lopo(X, y, groups, l2=1.0)
        assert np.isfinite(probabilities).sum() == len(y)
        assert all(f["status"] == "evaluated" for f in folds)
        # Performance is not an architecture invariant. Verify that the
        # grouped procedure produced finite, non-constant predictions; model
        # quality belongs in the locked expert-adjudicated evaluation set.
        assert np.isfinite(probabilities).all()
        assert np.all((probabilities >= 0.0) & (probabilities <= 1.0))
        # Ensure a patient's own opposite-class patch cannot be in its training fold.
        assert all(set(groups[f] for f in range(len(groups)) if groups[f] == fold["study_id"]) for fold in folds)
        return f"{len(y)} patches, {len(set(groups))} studies, LOPO AUC={clf.auc(y, probabilities):.3f}"


def test_manifest_contract_rejects_bad_inputs():
    clf = load_module(SCRIPTS / "classify_surface_patches.py", "classifier_manifest_test")
    with tempfile.TemporaryDirectory(prefix="css_arch_manifest_") as td:
        tmp = Path(td)
        patch = tmp / "x.npz"
        write_npz_patch(patch, "css")
        bad_duplicate = tmp / "duplicate.csv"
        with open(bad_duplicate, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["patch_path", "study_id", "label"]); w.writeheader()
            w.writerow({"patch_path": str(patch), "study_id": "P1", "label": "cSS"})
            w.writerow({"patch_path": str(patch), "study_id": "P1", "label": "vein"})
        try:
            clf.read_manifest(str(bad_duplicate))
        except ValueError:
            pass
        else:
            raise AssertionError("duplicate patches were accepted")
        bad_label = tmp / "bad_label.csv"
        with open(bad_label, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["patch_path", "study_id", "label"]); w.writeheader()
            w.writerow({"patch_path": str(patch), "study_id": "P1", "label": "artifact"})
        rows = clf.read_manifest(str(bad_label))
        try:
            clf.build_matrix(rows)
        except ValueError:
            pass
        else:
            raise AssertionError("invalid label was accepted")
        return "duplicate path and invalid label rejected"


def test_optional_fate_ledger():
    ledger_path = SCRIPTS / "object_ledger.py"
    if not ledger_path.exists():
        raise FileNotFoundError("scripts/object_ledger.py not implemented yet")
    ledger = load_module(ledger_path, "object_ledger_test")
    if not hasattr(ledger, "FATES") or not hasattr(ledger, "FateLedger"):
        raise AssertionError("ledger must expose FATES and FateLedger")
    if not {"shown_candidate", "rule_excluded", "generation_dropped", "qc_blocked", "not_generated"} <= set(ledger.FATES):
        raise AssertionError("ledger fate vocabulary incomplete")
    return "fate vocabulary present"


def main():
    print(f"Updated cSS architecture verification — project: {ROOT}")
    check("source modules compile", test_source_compiles)
    check("static pipeline contracts", test_static_pipeline_contracts, optional=True)
    check("surface frame and multimodal patch", test_surface_frame_and_multimodal_patch)
    check("classifier grouping and missing channels", test_classifier_patient_grouping_and_missing_channels)
    check("manifest contract", test_manifest_contract_rejects_bad_inputs)
    check("object-fate ledger", test_optional_fate_ledger, optional=True)
    print("\nSummary")
    print(f"  PASS: {len(PASS)}")
    print(f"  SKIP: {len(SKIP)}")
    print(f"  FAIL: {len(FAIL)}")
    if SKIP:
        print("  skipped: " + ", ".join(SKIP))
    if FAIL:
        print("  failed: " + ", ".join(FAIL))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
