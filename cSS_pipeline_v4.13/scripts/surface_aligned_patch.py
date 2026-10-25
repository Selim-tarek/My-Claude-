#!/usr/bin/env python3
"""Surface-aligned multimodal patches per detector object (v5.2, master plan Stage 2). Research software.
usage: surface_aligned_patch.py SUBJECT [--lists candidates,excluded] [--ids 3,7] [--size-mm 24,12,8]
                                [--spacing-mm 0.5,0.5,0.5]
For every object of review/ID_candidates|excluded.nii.gz a patch of fixed PHYSICAL size is sampled in a right-handed
frame at the object:
  u  the object's long axis (PCA, world mm) projected into the cortical tangent plane
  v  cross(n, u): across the surface
  n  outward pial normal (from tissue toward CSF): gradient of work/ID_pial_sd.nii.gz (detect_css.py), mapped to world
     coordinates with the inverse transpose of the affine (oblique data included), smoothed over 1 mm. The sign is
     checked against the labels (tissue must lie on the -n side); a reversed map is flipped and recorded.
Channels (each sampled on the same world grid, linear interpolation, NaN outside the image):
  swi    data/ID_swi.nii                         required
  phase  data/ID_phase.nii.gz                    optional; used only if on the SWI grid and newer than the SWI
  flair  work/ID_flair_swispace.nii.gz           optional; newer than the SWI
  labels work/ID_seg_swispace.nii.gz             optional; nearest neighbour; newer than the SWI
A missing or stale channel is NOT zero-filled: it is absent from the file, valid_<ch> is all False and the reason is
in the metadata. Intensity channels are valid only inside the image AND inside the brain mask (SWI > 0), so outside-
brain zeros never enter patch features. Objects with < 3 voxels are skipped with a reason.
Writes review/ID_patches/<list>_<id>.npz (image_<ch>, valid, valid_<ch>, world_grid_mm, center_world_mm, center_ijk,
u_world, v_world, n_world, metadata_json) and review/ID_patches.csv (one row per object, also the skipped ones).
Library entry: extract_patch(images, pial, mask, ...) on nibabel images (each channel sampled through its OWN affine;
the CLI only passes channels already on the SWI grid)."""
import sys, os, json, argparse, numpy as np, nibabel as nib, pandas as pd
from scipy import ndimage as ndi
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from css_common import base as _base, cand_uid, stale_against, build_id

PATCH_SCHEMA = "css-patch/1.0"


def unit(x, eps=1e-8):
    n = float(np.linalg.norm(x))
    if n < eps: raise ValueError("near-zero vector")
    return x / n


def to_world(aff, ijk):
    x = np.asarray(ijk, float).reshape(-1, 3)
    return (aff @ np.c_[x, np.ones(len(x))].T).T[:, :3].reshape(np.asarray(ijk).shape)


def to_voxel(aff, xyz):
    x = np.asarray(xyz, float).reshape(-1, 3)
    return (np.linalg.inv(aff) @ np.c_[x, np.ones(len(x))].T).T[:, :3].reshape(np.asarray(xyz).shape)


def world_normal(pial_sd, aff, ijk, vox, half_mm=6.0):
    """outward unit normal at voxel ijk: gradient of the (1 mm smoothed) signed distance in a local crop, as a world
    covector (inverse transpose of the affine's linear part)"""
    r = np.ceil(half_mm / np.asarray(vox)).astype(int); c = np.round(ijk).astype(int)
    lo = np.maximum(c - r, 0); hi = np.minimum(c + r + 1, pial_sd.shape)
    crop = ndi.gaussian_filter(pial_sd[tuple(slice(a, b) for a, b in zip(lo, hi))].astype(np.float64), 1.0 / np.asarray(vox))
    g = np.array([np.gradient(crop, axis=a)[tuple(c - lo)] for a in range(3)])        # per voxel index
    return unit(np.linalg.inv(np.asarray(aff, float)[:3, :3]).T @ g)


def frame(mask, pial_sd, aff, vox, tissue=None):
    ijk = np.argwhere(mask)
    if len(ijk) < 3: raise ValueError(f"object has {len(ijk)} voxels (< 3)")
    c_ijk = ijk.mean(0); near = ijk[np.argmin(np.linalg.norm(ijk - c_ijk, axis=1))]   # a voxel OF the object
    n = world_normal(pial_sd, aff, near, vox); flipped = False
    if tissue is not None:                                           # the labels decide which side is tissue
        side = lambda t: ndi.map_coordinates(tissue.astype(np.float32), to_voxel(aff, to_world(aff, near) + t * n).reshape(3, 1), order=0)[0]
        if side(+1.5) > 0.5 and side(-1.5) < 0.5: n, flipped = -n, True
    P = to_world(aff, ijk); C = P - P.mean(0)
    u = np.linalg.eigh(C.T @ C / max(len(P) - 1, 1))[1][:, -1]
    u = u - (u @ n) * n
    if np.linalg.norm(u) < 1e-5:
        e = np.eye(3)[np.argmin(np.abs(n))]; u = e - (e @ n) * n
    u = unit(u)
    if u[np.argmax(np.abs(u))] < 0: u = -u
    v = unit(np.cross(n, u)); u = unit(np.cross(v, n))
    return dict(center_ijk=c_ijk, center_world=to_world(aff, c_ijk), u=u, v=v, n=n, normal_flipped=flipped,
                pial_sd_at_center_mm=float(pial_sd[tuple(np.round(c_ijk).astype(int))]))


def grid(fr, size_mm, spacing_mm):
    ax = [(np.arange(int(round(s / d)) + 1) - int(round(s / d)) / 2) * d for s, d in zip(size_mm, spacing_mm)]
    U, V, N = np.meshgrid(*ax, indexing="ij")
    return fr["center_world"] + U[..., None] * fr["u"] + V[..., None] * fr["v"] + N[..., None] * fr["n"]


def sample(img_data, aff, xyz, order=1):
    co = np.moveaxis(to_voxel(aff, xyz), -1, 0)
    vals = ndi.map_coordinates(img_data, co, order=order, mode="constant", cval=np.nan, prefilter=False).astype(np.float32)
    inside = ndi.map_coordinates(np.ones(img_data.shape, np.float32), co, order=0, mode="constant", cval=0) > 0.5
    return vals, inside


_ARR = {}


def _arr(img):
    """float32 array of a nibabel image, loaded once per image object"""
    if id(img) not in _ARR: _ARR[id(img)] = (img, np.asarray(img.dataobj, np.float32))
    return _ARR[id(img)][1]


def extract_patch(images, pial, mask, size_mm=(24, 12, 8), spacing_mm=(0.5, 0.5, 0.5), obj_id=None, labels=None,
                  brain=None):
    """one surface-aligned patch.
    images: {channel: nibabel image}; pial: signed pial distance image (mm, + outward); mask: object label image
    (obj_id selects one object, None = every non-zero voxel); labels: optional segmentation image (tissue side check);
    brain: optional brain mask on the SWI grid (default: filled SWI > 0). Intensity channels are valid only inside
    their image and inside the brain; 'labels' only inside its image.
    -> dict(shape, channels, arrays {image_<ch>}, valids {valid_<ch>}, valid (= valid_swi, else any channel), frame, grid)
    raises ValueError for an object of < 3 voxels"""
    aff = pial.affine; vox = np.array(pial.header.get_zooms()[:3], float)
    M = _arr(mask); m = (M == obj_id) if obj_id is not None else (M > 0)
    tissue = None
    if labels is not None: L = _arr(labels); tissue = (L > 0) & (L != 24)
    fr = frame(m, _arr(pial), aff, vox, tissue)
    G = grid(fr, size_mm, spacing_mm)
    brain_s = None
    if "swi" in images:
        if brain is None: brain = ndi.binary_fill_holes(_arr(images["swi"]) > 0)
        brain_s = sample(brain.astype(np.float32), images["swi"].affine, G, order=0)[0] > 0.5
    arrays, valids = {}, {}
    for name, img in images.items():
        vals, inside = sample(_arr(img), img.affine, G, 0 if name == "labels" else 1)
        ok = inside & np.isfinite(vals)
        if name != "labels" and brain_s is not None: ok &= brain_s
        arrays[f"image_{name}"] = vals; valids[f"valid_{name}"] = ok
    valid = valids["valid_swi"] if "valid_swi" in valids else np.any(list(valids.values()), axis=0)
    return dict(shape=tuple(G.shape[:3]), channels=sorted(images), arrays=arrays, valids=valids, valid=valid, grid=G,
                frame=dict(u_world=fr["u"], v_world=fr["v"], n_world=fr["n"], center_world=fr["center_world"],
                           center_ijk=fr["center_ijk"], normal_flipped=fr["normal_flipped"],
                           pial_sd_at_center_mm=fr["pial_sd_at_center_mm"]))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("subject"); ap.add_argument("--lists", default="candidates,excluded"); ap.add_argument("--ids", default="")
    ap.add_argument("--size-mm", default="24,12,8"); ap.add_argument("--spacing-mm", default="0.5,0.5,0.5")
    a = ap.parse_args(); s = a.subject; B = _base()
    size = [float(x) for x in a.size_mm.split(",")]; spc = [float(x) for x in a.spacing_mm.split(",")]
    swi_p, psd_p = f"{B}/data/{s}_swi.nii", f"{B}/work/{s}_pial_sd.nii.gz"
    if not os.path.exists(psd_p): sys.exit(f"{psd_p} missing - re-run detect_css.py {s} (v5.2 writes it)")
    if stale_against(psd_p, [swi_p]): sys.exit(f"{psd_p} is older than the SWI - re-run detect_css.py {s}")
    swi = nib.load(swi_p); aff = swi.affine
    brain = ndi.binary_fill_holes(_arr(swi) > 0); psd_img = nib.load(psd_p)
    chans, missing = {"swi": swi}, {}
    for name, p in (("phase", f"{B}/data/{s}_phase.nii.gz"), ("flair", f"{B}/work/{s}_flair_swispace.nii.gz"),
                    ("labels", f"{B}/work/{s}_seg_swispace.nii.gz")):
        if not os.path.exists(p): missing[name] = "absent"; continue
        im = nib.load(p)
        if im.shape[:3] != swi.shape[:3] or not np.allclose(im.affine, aff, atol=1e-3): missing[name] = "not on the SWI grid"; continue
        if stale_against(p, [swi_p]): missing[name] = "older than the SWI"; continue
        chans[name] = im
    out_dir = f"{B}/review/{s}_patches"; os.makedirs(out_dir, exist_ok=True)
    ids = {int(x) for x in a.ids.split(",") if x.strip()}
    rows = []
    for lst in [x for x in a.lists.split(",") if x]:
        mp, cp = f"{B}/review/{s}_{lst}.nii.gz", f"{B}/review/{s}_{lst}.csv"
        if not (os.path.exists(mp) and os.path.exists(cp)): continue
        M_img = nib.load(mp); d = pd.read_csv(cp)
        idc = "cand_id" if lst == "candidates" else "excl_id"
        for _, r in d.iterrows():
            oid = int(r[idc])
            if ids and oid not in ids: continue
            row = dict(study_id=s, list=lst, obj_id=oid, cand_uid=cand_uid(s, r), patch_path="", status="skipped", skip_reason="")
            try:
                P = extract_patch(chans, psd_img, M_img, size, spc, obj_id=oid, labels=chans.get("labels"), brain=brain)
            except ValueError as e:
                row["skip_reason"] = str(e); rows.append(row); continue
            fr, G, arrays, valids, valid = P["frame"], P["grid"], P["arrays"], P["valids"], P["valid"]
            for name in missing: valids[f"valid_{name}"] = np.zeros(G.shape[:3], bool)
            meta = dict(schema=PATCH_SCHEMA, subject=s, list=lst, obj_id=oid, cand_uid=row["cand_uid"], size_mm=size,
                        spacing_mm=spc, axis_semantics={"0": "u_along_candidate", "1": "v_across_surface", "2": "n_outward"},
                        normal_sign="positive_outward_toward_csf", normal_flipped=fr["normal_flipped"],
                        pial_sd_at_center_mm=fr["pial_sd_at_center_mm"], channels=sorted(chans), missing_channels=missing,
                        build=build_id(), research_only=True)
            path = f"{out_dir}/{lst}_{oid}.npz"
            np.savez_compressed(path, **arrays, **valids, valid=valid, world_grid_mm=G.astype(np.float32),
                                center_world_mm=fr["center_world"], center_ijk=fr["center_ijk"], u_world=fr["u_world"],
                                v_world=fr["v_world"], n_world=fr["n_world"], metadata_json=json.dumps(meta))
            row.update(patch_path=os.path.relpath(path, B), status="ok", channels=";".join(sorted(chans)),
                       missing_channels=";".join(f"{k}:{v}" for k, v in missing.items()),
                       valid_frac_swi=round(float(valid.mean()), 3), normal_flipped=int(fr["normal_flipped"]),
                       pial_sd_at_center_mm=round(fr["pial_sd_at_center_mm"], 2))
            rows.append(row)
    pd.DataFrame(rows).to_csv(f"{B}/review/{s}_patches.csv", index=False)
    ok = sum(r["status"] == "ok" for r in rows)
    print(f"{s}: {ok} surface-aligned patches ({size} mm, spacing {spc} mm) -> {out_dir}; {len(rows) - ok} skipped; "
          f"channels {sorted(chans)}" + (f"; missing {missing}" if missing else "") + "  [research software]")
