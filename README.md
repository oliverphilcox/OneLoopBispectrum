# OneLoopBispectrum
Computation of the one-loop bispectrum of galaxies in redshift space, using Mathematica. This is described in [Philcox+22](COMING SOON). The code contains a number of parts:
1. *Mathematica:* Computation and simplification of the bispectrum kernels [cosmology-independent]
2. *Mathematica:* Computation of the bispectrum templates, given the simplified kernels 
3. *Mathematica:* Read-in and combination of the bispectrum templates, for efficient application
4. *Python:* Interpolation and extraction of the one-loop bispectrum
Steps 1 and 2 are performed for both flattened (k1 = k2 + k3) and full bispectrum shapes, to allow for efficient and numerically stable interpolation. Note that step 1 can be omitted unless the bispectrum kernels themselves are being modified.

The computed spectra can be interfaced with the tree-level bispectrum likelihoods found [here](https://github.com/oliverphilcox/full_shape_likelihoods). A necessary input to stage (2) is the IR-resummed power spectrum; this can be computed using [CLASS-PT](https://github.com/michalychforever/CLASS-PT).

### Authors
- [Oliver Philcox](mailto:ohep2@cantab.ac.uk) (Princeton / IAS) 
