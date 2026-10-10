"""Labels on the SWI grid -> work/ID_seg_swispace.nii.gz (masked to SWI coverage).
Source, best first:
  work/ID_t1seg_swispace.nii.gz  T1-derived labels already in SWI space (prep_anat.sh, v4)
  synthseg/ID_seg.nii.gz         SynthSeg run on the SWI itself (v3 behaviour)"""
import signal; signal.signal(signal.SIGPIPE, signal.SIG_DFL)   # quiet when piped to head
import sys, os, numpy as np, nibabel as nib
from nibabel.processing import resample_from_to
from scipy import ndimage as ndi
subj = sys.argv[1]
base = os.environ.get("CSS_BASE", os.path.expanduser("~/css_project"))
swi = nib.load(f"{base}/data/{subj}_swi.nii")
t1p = f"{base}/work/{subj}_t1seg_swispace.nii.gz"
src = t1p if os.path.exists(t1p) and os.path.getmtime(t1p) >= os.path.getmtime(f"{base}/data/{subj}_swi.nii") \
    else f"{base}/synthseg/{subj}_seg.nii.gz"
seg = nib.load(src)
seg_swi = resample_from_to(seg, swi, order=0).get_fdata().astype(np.int32)
n_before = (seg_swi > 0).sum()
seg_swi[~ndi.binary_fill_holes(swi.get_fdata() > 0)] = 0
out = nib.Nifti1Image(seg_swi, swi.affine, swi.header); out.set_data_dtype(np.int32)
nib.save(out, f"{base}/work/{subj}_seg_swispace.nii.gz")
print(f"{subj}: labels from {os.path.basename(src)}; labelled voxels before mask {n_before:,}, "
      f"after mask {(seg_swi>0).sum():,}")
