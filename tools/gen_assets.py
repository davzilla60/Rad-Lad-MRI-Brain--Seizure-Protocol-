"""Build ../assets for the MRI EPILEPSY (hippocampus) simulator.

    python -I tools/gen_assets.py "<.../ds005602/sub-2/anat>"

Source: IDEAS - Imaging Database for Epilepsy And Surgery (OpenNeuro ds005602, CC0; Taylor et al. 2024), sub-2:
right hippocampal sclerosis (ILAE type 1 on histology), GE Discovery MR750 3T, UCLH. 3D FLAIR (1 mm isotropic) and
3D T1 (0.94 x 0.94 x 1.5 mm). Working frame = the FLAIR's scanner space (LPS mm). The T1 is moved into it by the rigid
registration in tools/t1_to_flair.json (made by epilepsy-mri/tools/register_t1.py: ~3.7 mm, ~2 deg).
  * SAG T1: reformatted from the 3D T1, parallel to the falx.
  * AX FLAIR: reformatted from the 3D FLAIR, parallel to the long axis of the hippocampus.
  * COR FLAIR, COR T1: reformatted perpendicular to the long axis of the hippocampus.
  * LOCALISER: synthesised from the 3D T1 (8 mm slab, 256 matrix, 280 mm FOV, coil shading, noise).
"""
import sys, os, glob, json
import numpy as np, nibabel as nib, cv2
from scipy.ndimage import map_coordinates, gaussian_filter1d

ANAT = sys.argv[1]
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "assets"); QA = os.path.join(HERE, "..", "qa")
os.makedirs(OUT, exist_ok=True); os.makedirs(QA, exist_ok=True)
rng = np.random.default_rng(11)
def unit(a): a = np.asarray(a, float); return a / np.linalg.norm(a)
FLIP = np.array([-1.0, -1.0, 1.0])

class Vol:
    def __init__(self, f, pre=None):
        img = nib.load(f); self.v = np.asanyarray(img.dataobj).astype(np.float32)
        aff = img.affine if pre is None else pre @ img.affine
        self.inv = np.linalg.inv(aff)
    def sample(self, P):
        ijk = (P * FLIP) @ self.inv[:3, :3].T + self.inv[:3, 3]
        return map_coordinates(self.v, [ijk[..., 0], ijk[..., 1], ijk[..., 2]], order=1, mode="constant", cval=0.0)
M = np.array(json.load(open(os.path.join(HERE, "t1_to_flair.json")))["M_t1_to_flair_ras"])
FL = Vol(glob.glob(os.path.join(ANAT, "*_FLAIR.nii.gz"))[0])
T1 = Vol(glob.glob(os.path.join(ANAT, "*_T1w.nii.gz"))[0], M)
HI = {"fl": float(np.percentile(FL.v[FL.v > np.percentile(FL.v, 40)], 99.7)),
      "t1": float(np.percentile(T1.v[T1.v > np.percentile(T1.v, 40)], 99.7))}

# ------------------------------------------------------------------ anatomy (LPS mm, picked on qa renders, ±3 mm)
HIP_ANG = -18.0                    # hippocampal long axis on the sagittal: rises ~18 deg toward the back (both sides)
a = np.radians(HIP_ANG)
VH = np.array([0.0, np.cos(a), -np.sin(a)])         # along the hippocampus (+posterior, rising)
UA = np.array([1.0, 0, 0]); NA = np.cross(UA, VH)   # axial parallel to the hippocampus
UC, VC, NC = UA, -NA, VH                            # coronal perpendicular to it
US, VS, NS = np.array([0, 1.0, 0]), np.array([0, 0, -1.0]), np.array([-1.0, 0, 0])
BC = np.array([3.0, -25.0, 0.0])
print("axial normal", NA.round(3), " coronal normal", NC.round(3))

def to8(x, hi): return np.clip(x / hi * 255, 0, 255).astype(np.uint8)
calib = {"ser": {}}
def reformat(sid, vol, key, plane, label, params, u, v, n, ctr, N, fov, thk, sp, count, spl):
    ps = fov / N
    I, J = np.meshgrid(np.arange(N) - (N - 1) / 2, np.arange(N) - (N - 1) / 2)
    first = ctr - n * sp * (count - 1) / 2; ipps = []
    for k in range(count):
        c = first + n * sp * k
        base = c + I[..., None] * ps * u + J[..., None] * ps * v
        acc = np.mean([vol.sample(base + n * t) for t in np.linspace(-thk / 2 + thk / 10, thk / 2 - thk / 10, 4)], axis=0)
        cv2.imwrite(os.path.join(OUT, "%s_%02d.jpg" % (sid, k)), to8(acc, HI[key]), [cv2.IMWRITE_JPEG_QUALITY, 80])
        ipps.append(c - u * ps * (N - 1) / 2 - v * ps * (N - 1) / 2)
    calib["ser"][sid] = {"task": plane, "plane": plane, "seq": label, "params": params, "n": count,
                         "ipp": [[round(float(x), 2) for x in p] for p in ipps], "row": [round(float(x), 5) for x in u],
                         "col": [round(float(x), 5) for x in v], "ps": round(ps, 6), "cols": N, "rows": N,
                         "nrm": [round(float(x), 5) for x in n], "sp": sp, "thk": thk, "spl": spl, "pad": 2}
    print(sid, count, "slices")

reformat("sag_t1", T1, "t1", "sag", "T1", "3D T1 · reformatted 3 mm", US, VS, NS, BC.copy(), 384, 240, 3.0, 3.0, 54, 3.0)
reformat("ax_flair", FL, "fl", "ax", "FLAIR", "3D FLAIR · reformatted 3 mm, along the hippocampus", UA, VH, NA, BC.copy(), 384, 240, 3.0, 3.0, 52, 4.0)
reformat("cor_flair", FL, "fl", "cor", "FLAIR", "3D FLAIR · reformatted 3 mm, perpendicular to the hippocampus", UC, VC, NC, BC + NC * (-8.0), 384, 240, 3.0, 3.0, 62, 4.0)
reformat("cor_t1", T1, "t1", "cor", "T1", "3D T1 · reformatted 3 mm, perpendicular to the hippocampus", UC, VC, NC, BC + NC * (-8.0), 384, 240, 3.0, 3.0, 62, 3.0)

# ------------------------------------------------------------------ synthetic localiser from the registered T1
LOC_FOV, LN = 280.0, 256; LP = LOC_FOV / LN
LOC_C = np.array([3.0, -25.0, -5.0])
LOC = {"sag": (US, VS, np.array([0.0, 0, 0])), "cor": (UA, np.array([0, 0, -1.0]), np.array([0, 0.0, 0])),
       "ax": (UA, np.array([0, 1.0, 0]), np.array([0, 0, 10.0]))}
calib["loc"] = {"n": LN, "fov": LOC_FOV}
for p, (r, c, off) in LOC.items():
    n = np.cross(r, c); ctr = LOC_C + off; SS = 2; N2 = LN * SS
    I, J = np.meshgrid((np.arange(N2) + 0.5) / SS - LN / 2, (np.arange(N2) + 0.5) / SS - LN / 2)
    base = ctr + I[..., None] * r * LP + J[..., None] * c * LP
    im = np.mean([T1.sample(base + n * t) for t in np.linspace(-3.5, 3.5, 6)], axis=0)
    im = cv2.resize(im.astype(np.float32), (LN, LN), interpolation=cv2.INTER_AREA)
    Pc = ctr + ((np.arange(LN) + 0.5) - LN / 2)[None, :, None] * r * LP + ((np.arange(LN) + 0.5) - LN / 2)[:, None, None] * c * LP
    rr = np.linalg.norm((Pc - BC) / [80, 100, 95], axis=-1)
    im = im * (0.8 + 0.3 * np.clip(rr, 0, 1.2) ** 2)
    im = gaussian_filter1d(im, 0.8, axis=0); sig = 0.025 * np.percentile(im, 99.5)
    im = np.hypot(im + rng.normal(0, sig, im.shape), rng.normal(0, sig, im.shape))
    hi = np.percentile(im[im > 0], 99.6)
    cv2.imwrite(os.path.join(OUT, "loc_" + p + ".png"), np.clip((im / hi) ** 0.9 * 255, 0, 255).astype(np.uint8))
    ipp = ctr - r * LP * LN / 2 - c * LP * LN / 2 + (r + c) * LP / 2
    calib["loc"][p] = {"ipp": [round(float(x), 2) for x in ipp], "row": r.tolist(), "col": c.tolist(), "ps": round(LP, 6), "cols": LN, "rows": LN}

calib["lm"] = {
    "conf":   BC.tolist(),
    "rhh":    [-28.0, -33.0, -27.0],   # right hippocampal head
    "rht":    [-28.0, -9.0, -19.0],    # right hippocampal tail
    "rhb":    [-28.0, -21.0, -23.0],   # right hippocampal body
    "lhh":    [28.0, -35.0, -27.0],    # left hippocampal head
    "lht":    [28.0, -7.0, -17.0],     # left hippocampal tail
    "lhb":    [28.0, -21.0, -22.0],    # left hippocampal body
    "rtp":    [-35.0, -58.0, -28.0],   # right temporal pole
    "ltp":    [35.0, -60.0, -28.0],    # left temporal pole
    "vertex": [3.0, -30.0, 48.0],
    "fm":     [3.0, -18.0, -50.0],
    "cbl":    [3.0, 15.0, -43.0],
    "fpole":  [3.0, -103.0, 2.0],
    "opole":  [3.0, 53.0, -8.0],
    "rtl":    [-58.0, -20.0, -25.0],   # right temporal lobe, lateral
    "ltl":    [66.0, -20.0, -25.0],    # left temporal lobe, lateral
    "mid":    [3.0, -25.0, 10.0],
    "front":  [3.0, -118.0, 5.0], "back": [3.0, 68.0, -5.0],
    "rskull": [-74.0, -20.0, -5.0], "lskull": [84.0, -20.0, -5.0],
}
calib["real"] = {"ax": HIP_ANG, "cor": HIP_ANG + 90, "sag": 90.0}
json.dump(calib, open(os.path.join(OUT, "calib.json"), "w"), indent=1)
open(os.path.join(OUT, "calib.js"), "w").write("window.MRS_CALIB=" + json.dumps(calib, separators=(",", ":")) + ";\n")

def qa(sid, keys):
    s = calib["ser"][sid]; ipp = np.array(s["ipp"]); r = np.array(s["row"]); c = np.array(s["col"]); nrm = np.array(s["nrm"])
    tiles = []
    for key in keys:
        P = np.array(calib["lm"][key]); k = int(np.argmin(np.abs((P - ipp) @ nrm)))
        im = cv2.cvtColor(cv2.imread(os.path.join(OUT, "%s_%02d.jpg" % (sid, k)), 0), cv2.COLOR_GRAY2BGR)
        x, y = (P - ipp[k]) @ r / s["ps"], (P - ipp[k]) @ c / s["ps"]
        cv2.circle(im, (int(x), int(y)), 7, (0, 0, 255), 1); cv2.putText(im, key, (int(x) + 9, int(y)), 0, 0.45, (0, 255, 255), 1)
        tiles.append(im)
    cv2.imwrite(os.path.join(QA, "qa_" + sid + ".png"), np.hstack(tiles))
qa("sag_t1", ["rhb", "lhb", "vertex", "fm", "fpole", "opole"])
qa("ax_flair", ["rhh", "lhb", "rtp", "vertex", "cbl", "front"])
qa("cor_flair", ["rhb", "rhh", "rht", "rtp", "fpole", "opole"])
qa("cor_t1", ["rhb", "rhh", "rht", "rtp"])
print("ok")
