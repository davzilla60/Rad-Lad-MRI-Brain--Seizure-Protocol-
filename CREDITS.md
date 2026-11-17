# Credits

## Imaging

Adapted from **IDEAS — the Imaging Database for Epilepsy And Surgery**, OpenNeuro
[ds005602](https://openneuro.org/datasets/ds005602), released under **CC0 1.0** (public domain).

> Taylor, P. N., Wang, Y., Simpson, C., Janiukstyte, V., Horsley, J., Leiberg, K., Little, B., Clifford, H.,
> Adler, S., Vos, S. B., Winston, G. P., McEvoy, A. W., Miserocchi, A., de Tisi, J., & Duncan, J. S. (2024).
> The Imaging Database for Epilepsy And Surgery (IDEAS). https://arxiv.org/abs/2406.06731

Participant **sub-2**: scanned at UCLH / the National Hospital for Neurology and Neurosurgery, GE Discovery MR750 3 T.
The pathology (right hippocampal sclerosis, ILAE type 1) and the outcome (seizure-free) come from the IDEAS clinical
metadata table. The wording of the request in the simulator is a paraphrase, not the original referral.

### Changes made

- The 3D T1 was rigidly registered to the 3D FLAIR (≈ 3.7 mm, 2°) so both share one space.
- Sagittal T1, axial FLAIR and coronal FLAIR / T1 are **reformatted** (3 mm slabs) from the 3D volumes: sagittal
  parallel to the falx, axial parallel to the long axis of the hippocampus, coronal perpendicular to it.
- The 3-plane localiser is **synthesised** from the 3D T1 (8 mm slab, 280 mm FOV, coil shading, noise).
- Pixel values windowed to 8-bit JPEG / PNG.
- The 3-D room shows a generic closed-bore magnet and head coil, not any real product.

Educational simulation — not a medical device, not for clinical use.
