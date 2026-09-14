# Recovery of the FLAMINGO Stacked kSZ Profile

The aim of this repository is to reproduce the stacked kSZ profile from the public FLAMINGO data. This will form the validation data vector against which to build the forward model. The measurement follows McCarthy et al. 2024 (arXiv:2410.19905).

## Requirements

```bash
python -m pip install -r requirements.txt
```

## Quick Start

```bash
python 00_synthetic_test.py  # verify the machinery, no data needed
python 01_explore_flamingo.py  # find out what is in the files
python 02_build_ksz_map.py --shells 9 10 11 12 --remote
python 03_select_sample.py --sample cmass --halo_dir <dir with lightcone_halos_*.hdf5>
python 04_stack_cap.py --map dT_ksz_*.npy --catalogue sample_cmass_*.npz --beam 1.3
python 05_plot_datavector.py --cap cap_profile.npz --label "FLAMINGO L1_m9"
```

## Pipeline

### ksz_utils.py:
- Contains the physics and geometry that the other scripts import.

### 00_synethetic_test.py:
- Constructs an artificial sky with a signal of known amplitude and scape dependence, and checks that the framework recovers it.
- The sanity checks performed are:
  - Indexing logic of the CAP profile;
  - Normalisation of the CAP profile;
  - Integration to zero of a constant sky;
  - Recovery of the injected signal from the velocity weighted stack.

### 01_explore_flamingo.py:
- Reads the FLAMINGO metadata from the Durham server (Cosma) and prints:
  - Which simulations exist and what their directories are called;
  - The layout of the shell files and which maps are inside each one;
  - The redshift and comoving bounds of each shell, and the pixel size;
  - The halo-lightcone file list and redshift coverage;
  - The SOAP dataset names.

### 02_build_ksz_map.py:
- Builds the kSZ temperature map from the Doppler B parameter, $b$, using $\Delta T = −T_{\rm CMB} \times \Sigma_{i \in {\rm shells}} b_i$.

### 03_select_sample.py:
- Build the sample to stack on. There are two selection modes:
  - Galaxy-like. This follows McCarthy et al. 2024 in which a minimum stellar mass cut is applied.
    ```bash
    python 03_select_sample.py --sample cmass
    ```
  - Halo-mass bin. This selects haloes directly in a bin of $M_{500{\rm c}}$ and stacks them.
    ```bash
    python 03_select_sample.py --halo-mass-bin 13.0 13.5
    ```

### 04_stack_cap.py:
- Computes the final stacked kSZ measurement via CAP on each galaxy, with velocity weighting, and obtains bootstrap errors.
