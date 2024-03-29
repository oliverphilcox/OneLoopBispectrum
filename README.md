# One Loop Galaxy Bispectrum
Computation of the one-loop bispectrum of galaxies in redshift space, using Mathematica. This is described in [Philcox+22](https://arxiv.org/abs/2206.02800). The code contains a number of components:
1. [Kernels.nb](Kernels.nb) \& [KernelsFlat.nb](KernelsFlat.nb): Computation and simplification of the bispectrum kernels using Mathematica [cosmology-independent]. These are saved in the [```kernels```](kernels/) directory.
2. [Spectra.wls](Spectra.wls) \& [SpectraFlat.wls](SpectraFlat.wls): Computation of the full and flattened bispectrum templates, given the pre-computed integration kernels. These scripts compute the templates for all bias parameters and a single angular component. They are designed to be run on a HPC cluster using Mathematica. We provide a sample [SLURM submission script](ComputeSpectra.slurm) for computing all templates.
3. [Output.wls](Output.wls): Read-in and combination of the bispectrum templates using Mathematica. These are saved as two HDF5 files.
4. [Bispectrum Interpolation.ipynb](Bispectrum%20Interpolation.ipynb): Interpolation and extraction of the one-loop bispectrum using Python. This demonstrates the practical usage of the one-loop bispectrum, including full radial and angular bin-integration, and Alcock-Paczynski distortions.
Steps 1-3 are performed for both flattened (k1 = k2 + k3) and full bispectrum shapes, to allow for efficient and numerically stable interpolation. Note that step 1 only needs to be run if the bispectrum kernels themselves are being modified, else the default kernels provided can be used.

A necessary input to stage (2) is the IR-resummed power spectrum in the [```pk```](pk/) directory; this can be computed using [CLASS-PT](https://github.com/michalychforever/CLASS-PT), and we provide a sample spectrum. High-resolution bispectrum templates (the output of step 3) computed using the PT Challenge best-fit power spectrum can be downloaded from [Google Drive](https://drive.google.com/drive/folders/1SVbW8oaMYxpCuLLlnd1dIcrg0vfp-wi_?usp=sharing), and are those used in the [Philcox+22](https://arxiv.org/abs/2206.02800) paper.

The computed spectra can be interfaced with the power spectrum and bispectrum likelihoods found [here](https://github.com/oliverphilcox/full_shape_likelihoods). 

### Authors
- [Oliver Philcox](mailto:ohep2@cantab.ac.uk) (Princeton / IAS) 
- Mikhail Ivanov (IAS)
- Marko Simonovic (CERN)
