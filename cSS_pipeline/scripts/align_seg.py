import signal; signal.signal(signal.SIGPIPE, signal.SIG_DFL)   # quiet when piped to head
import sys, os, numpy as np, nibabel as nib
from nibabel.processing import resample_from_to
from scipy import ndimage as ndi
subj = sys.argv[1]
base = os.environ.get("CSS_BASE", os.path.expanduser("~/css_project"))
swi = nib.load(f"{base}/data/{subj}_swi.nii")
seg = nib.load(f"{base}/synthseg/{subj}_seg.nii.gz")
seg_swi = resample_from_to(seg, swi, order=0).get_fdata().astype(np.int32)
n_before = (seg_swi > 0).sum()
seg_swi[~ndi.binary_fill_holes(swi.get_fdata() > 0)] = 0
out = nib.Nifti1Image(seg_swi, swi.affine, swi.header); out.set_data_dtype(np.int32)
nib.save(out, f"{base}/work/{subj}_seg_swispace.nii.gz")
print(f"{subj}: labelled voxels before mask {n_before:,}, after mask {(seg_swi>0).sum():,}")
