import argparse
import numpy as np
import healpy as hp

from ksz_utils import (cap_profiles_catalogue, velocity_weighted_stack,
                       bootstrap_errors, gaussian_beam_smooth,
                       SR_TO_ARCMIN2)

# --------------------------------------------------------------------------
# Stack the map around the sample with the CAP filter
# --------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--map", required=True)  # map from 02_build_map.py
    p.add_argument("--catalogue", required=True)  # catalogue from 03_select_sample.py 
    p.add_argument("--theta-min", type=float, default=1.0)
    p.add_argument("--theta-max", type=float, default=6.0)
    p.add_argument("--n-ap", type=int, default=8)
    p.add_argument("--theta-list", type=float, nargs="+", default=None,
                   help="explicit aperture radii in arcmin, overriding the "
                        "min/max/n grid. Schaan+21 (which McCarthy+24 follow) "
                        "use ~1 to 6 arcmin; pass their exact grid here when "
                        "comparing point by point")
    p.add_argument("--beam-file", default=None,
                   help="text file with the beam transfer function b_l, one "
                        "value per multipole starting at l=0 (e.g. the measured "
                        "ACT beam). Overrides the Gaussian --beam")
    p.add_argument("--beam", type=float, default=0.0,
                   help="Gaussian beam FWHM in arcmin (ACT f150 ~1.3, f090 ~2.1); "
                        "0 disables")
    p.add_argument("--max-objects", type=int, default=0,
                   help="subsample for a quick test run")
    p.add_argument("--n-boot", type=int, default=500)
    p.add_argument("--out", default="cap_profile.npz")
    args = p.parse_args()

    sky = np.load(args.map)
    npix = sky.size
    nside = int(round(np.sqrt(npix / 12)))
    print(f"map: nside {nside}, {npix} pixels, rms {sky.std():.3f} uK")

    # ----------------------------------------------------------------
    # Step 1: Convolve with beam
    # ----------------------------------------------------------------
    if args.beam_file:
        # McCarthy+24 use the measured ACT beam, applied in multipole space with
        # healpy's almxfl. That is the same operation as a Gaussian smooth with
        # a different b_l, so swapping one for the other is a two-line change.
        bl = np.loadtxt(args.beam_file)
        bl = bl[:, -1] if bl.ndim > 1 else bl
        lmax = min(len(bl) - 1, 3 * nside - 1)
        print(f"convolving with the beam in {args.beam_file}, lmax={lmax}")
        alm = hp.map2alm(sky.astype(np.float64), lmax=lmax)
        alm = hp.almxfl(alm, bl[:lmax + 1])
        sky = hp.alm2map(alm, nside=nside, lmax=lmax).astype(np.float32)
    elif args.beam > 0:
        print(f"convolving with a {args.beam}' FWHM Gaussian beam "
              "(this is the expensive step; lmax=3*nside-1)")
        sky = gaussian_beam_smooth(sky.astype(np.float64), args.beam).astype(np.float32)

    # ----------------------------------------------------------------
    # Step 2: Apply CAP filter
    # ----------------------------------------------------------------
    cat = np.load(args.catalogue)
    vec, v_r = cat["unit_vec"], cat["v_r_kms"]
    if args.max_objects and len(vec) > args.max_objects:
        idx = np.random.default_rng(0).choice(len(vec), args.max_objects, replace=False)
        vec, v_r = vec[idx], v_r[idx]
    print(f"catalogue: {len(vec)} objects, "
          f"satellite fraction {1 - cat['is_central'].mean():.3f}")

    theta = (np.asarray(args.theta_list, dtype=float) if args.theta_list
             else np.geomspace(args.theta_min, args.theta_max, args.n_ap))
    theta_rad = np.radians(theta / 60.0)
    print("apertures [arcmin]:", np.round(theta, 2))

    print("\nrunning CAP filter...", flush=True)
    T_cap = cap_profiles_catalogue(sky, nside, vec, theta_rad) * SR_TO_ARCMIN2

    # ----------------------------------------------------------------
    # Step 3: Stack and bootstrap
    # ----------------------------------------------------------------

    stack = velocity_weighted_stack(T_cap, v_r)
    _, err = bootstrap_errors(T_cap, v_r, n_boot=args.n_boot)

    # Null test: randomise the velocity signs. Should be consistent with zero.
    rng = np.random.default_rng(1)
    null = velocity_weighted_stack(T_cap, v_r * rng.choice([-1, 1], len(v_r)))

    print("\n theta[']    T_kSZ [uK arcmin^2]        null")
    for t, s, e, n in zip(theta, stack, err, null):
        print(f"  {t:6.2f}    {s:9.4f} +/- {e:6.4f}    {n:9.4f}")

    np.savez(args.out, theta_arcmin=theta, T_ksz=stack, T_err=err, T_null=null,
             n_obj=len(vec), beam_fwhm_arcmin=args.beam)
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
