# MR-I Ready — Epilepsy (hippocampus)

An embeddable, single-file teaching simulator for the Rad-Lad website, built on the MR-I Ready engine (forked from
the brain simulator). Every diagnostic image is reformatted from **one real 3 T MRI of a patient with surgically
proven right hippocampal sclerosis** (public-domain research data — see `CREDITS.md`) and is shown with its true
geometry. **It's a demonstration only.**

Marker: **Same Scan We Do Every Night** — sagittal T1 parallel to the falx with both hippocampi in the stack; axial parallel
to the long axis of the hippocampus and coronal perpendicular to it, both whole brain and in thin slices; and the
head straight, chin neutral, padded in the coil.

## Flow

1. **Patient set-up.** Supine, head first, 32-channel head coil. Choose the head position: straight and padded
   (correct), chin tipped up, or rolled to one side.
2. **Localiser.** A 3-plane localiser, synthesised from the patient's 3D T1, is "acquired" line by line from k-space.
3. **Sagittal T1 first.** Parallel to the falx, temporal lobe to temporal lobe (it must hold both hippocampi).
4. **Axial + coronal** unlock once the sagittal is scanned. Drag the sagittal T1 onto the **SAGITTAL** panel: it
   opens through the right hippocampus. Axial parallel to the hippocampal long axis (reference −18°: it rises toward
   the back, so a straight axial fails), coronal perpendicular (72°), graded ±8°.
5. **Scan → images.** Sagittal T1, axial FLAIR, coronal FLAIR + T1 (3 mm reformats).
6. **Debrief.** Coverage and angle per plan; hero images (sagittal with your lines against the hippocampal axis, and
   a coronal FLAIR zoomed on both hippocampi); the finding (right HS); a sequences table (HARNESS-MRI); the marker.

CLI: `neutral | chinup | tilted`, `load`, `loc`, `plan sag|ax|cor`, `auto`, `slices`, `thick`, `gap`, `fov`,
`angle`, `use <series>`, `restore <sag|cor|ax>`, `scan`, `back`, `finish`, `help`.

## Files

- `epilepsy.embed.html` is the self-contained deliverable; `index.html` is an identical copy.
- `epilepsy.html` + `app.js` + `room3d.js` + `assets/` are the editable sources. Rebuild with `python inline.py`
  (bump the `?v=` tags in `epilepsy.html` and `inline.py` together).
- `tools/gen_assets.py "<.../ds005602/sub-2/anat>"` rebuilds `assets/`. It uses `tools/t1_to_flair.json`, the rigid
  T1→FLAIR registration from `epilepsy-mri/tools/register_t1.py` (the two scans are ~3.7 mm / 2° apart).
