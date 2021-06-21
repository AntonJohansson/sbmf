# Slightly Beyond Mean-Field -- or something

Extends the Gross-Piatevskii (GP) equation using third order Rayleigh-Schrödinger or Epstein-Nesbet perturbation theory, for the more efficient study of Bose-Einstein condensate ground states. This is done via firstly solving the GP equation using the Self-Consistent Field (SCF) method in the Harmonic Oscillator (HO) basis, and then applying perturbation theory to the GP eigenstates. The SCF implementation works in an arbitrary basis whilst the perturbative calculations require are hard-coded to use the HO basis (a lot of integrals were precomputed in this basis).

## Building
```
$ make release
$ make debug
```
for release/debug mode, will produce ```build/sbmf.a```. Dependencies are currently handled in ```third_party``` and will be built with the project. Main dependencies are [OpenBlas](https://github.com/xianyi/openblas) for CBLAS/LAPACKE and [ArpackNG](https://github.com/opencollab/arpack-ng), these will be fetched and built from their respective git repositories using the ```third_party/build.sh``` script. CUDA is also required, and the system version will be used (assuming it is installed to ```usr/local/cuda```). TODO: more configuration options here.

## Reading the source

The code is structured in a unity-build manner, all public interface is available in ```include/sbmf/sbmf.h```, and ```src/sbmf/sbmf.c``` acts as the main source file.
