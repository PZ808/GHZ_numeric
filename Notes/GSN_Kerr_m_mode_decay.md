# SN circular Kerr seed: m-mode radial-decay check

Run on 2026-09-24. The asymptotic homogeneous Hertz completion plus its Bondi
residual gauge term has the expected **1/r decay in nn, nm, and mm** for the
circular prograde equatorial orbit below. This tests the asymptotic completion,
not the full sourced metric, a sum over all m, or the static sector.

## Parameters and data path

- M = mu = 1, a = 0.5, p = r0 = 10, e = 0, x = 1.
- m = +2 and -2, n = k = 0; omega_2 = 0.062261118482988327.
- GSN source modes ell = 2 through 8; spheroidal angular functions expanded
  into spherical coefficients solely to transfer the summed m-mode to C++.
- Julia GeneralizedSasakiNakamura v0.9.0, installed source previously checked
  against the local checkout. The exporter records the package path in the CSV.
- Example: Zinf_(2,2) = -0.0010444737265354943 + 0.0002838086433268029 i.
- Input: `Data/gsn_circular_a0.5_m2.csv`.
- Reduced fields use P_(m,s) = (1-z)^(|m+s|/2) (1+z)^(|m-s|/2).

The physical psi4 is rho^4 times the spin -2 Teukolsky variable. With the
paper's expansion in rho rather than 1/r, q = psi4^(1 circ) = -Zinf S.
This radial-expansion sign is separate from metric-signature conventions; the
metric formulas here follow the mostly-minus reconstruction and the existing
Schwarzschild normalization. See `GSN_Bondi_interface.md` for that interface.

The angular solve is genuinely on the LGL **m-mode grid**, using the existing
pole-factorized Held operators at nonzero a. No ell-by-ell seed inversion or
metric reconstruction is performed. We solve the equivalent coupled fourth-order
system rather than multiply it into a more ill-conditioned eighth-order system:

    A f + 3 i M omega fbar = (2 i / omega) q
    B fbar + 3 i M omega f = (2 i / omega) qbar
    A = Held-edth-prime^4, B = Held-edth^4.

Here qbar_m is the pointwise conjugate of q_-m, including its opposite frequency;
it is not the conjugate of q_m. The two reduced spin fields are independent
unknowns in this linear system. Row/column equilibration and iterative refinement
are used. For Nz=13 the relative matrix residual is 6.24e-11.

## Metric and the Kerr cancellation

`tools/bondi_kerr_m_mode_check.cpp` evaluates the full rational dependence on
rho = -1/(r-i a z), with analytic radial derivatives through order two and
spectral angular derivatives. It implements Sdagger (H.2b), the radial gauge
vector (E.3), and the Lie components (E.1d-f) from
[Hollands and Toomani](https://arxiv.org/abs/2403.20311).
The Hertz coefficients follow the triangular system in the local Bondi/GHZ
paper draft. Half of that paper's potential with Sdagger(Phi)+conjugate matches
the existing code's real-part normalization.

Use the explicit gauge construction **(J.1d)**:

    zeta_n^circ = Re[(Held-edth-prime Held-edth - rho_prime^circ) zeta_l^circ
                    - Omega^circ Held-edth-prime zeta_m^circ
                    + 2 bartau^circ zeta_m^circ].

In particular, the final tau term must be included. Omitting it in this
implementation leaves an O(1) remainder in h_nn for Kerr, although the a=0
control passes. In the paired m-mode expression its contribution is
bartau^circ zeta_m^circ + tau^circ zeta_mbar^circ. This is why testing only the
Schwarzschild formulas with a Kerr seed is insufficient.

The adjusted field is evaluated as reconstructed + Lie at each radius, without
manually deleting any growing or constant coefficient. All magnitudes in the
radial CSV and plot are maxima over the angular grid after restoring pole
factors. Complex equatorial samples are also saved for resolution comparisons.
The small resolution dependence of the sampled maximum need not equal the
convergence of the field at a common angular position.

## Results

Log-log secant slopes between r=200 and r=2000, Nz=13:

| Component | Kerr a=0.5 | Schwarzschild a=0 |
|---|---:|---:|
| nn | -1.00105158 | -1.00116929 |
| nm | -1.00497934 | -1.00505080 |
| mm | -0.99999866 | -1.00000000 |

The radius scan covers 20 to 20,000. At r=20,000, r times the maximum adjusted
nn amplitude is approximately 0.512513; it is 0.512514 at r=2000.

Validation:

- An independently generated a=0 SN seed was passed through this generic
  diagnostic. Both reconstruction and Lie pieces separately agree with the
  existing Schwarzschild Laurent implementation to 1.37e-13 relative error
  (r=20,100,500), before their cancellation.
- Nz=13,17,21: maximum relative differences of complex equatorial adjusted fields
  over the entire scan are below 3.8e-8 against Nz=13. The seed values at the
  equator agree to roughly 4e-13 relative. The coupled matrix residual grows
  to 6.34e-9 at Nz=21, so increasing Nz indefinitely is not beneficial in double
  precision.
- Separate GSN exports with source ell_max=6 and 8: maximum relative equatorial
  differences are 5.99e-7 (nn), 2.71e-6 (nm), 9.83e-8 (mm). These compare two
  cutoffs; they are not rigorous bounds on the omitted infinite tail.
- Four relevant CTest tests passed, including the new Kerr decay and generic
  Schwarzschild control tests, and the existing Schwarzschild data-collocation
  and metric-falloff tests. New tests require each slope to lie within 0.03 of
  -1 and the seed residual below 1e-7. The a=0 oracle tolerance is 1e-6.

## Reproduce

From the repository root, with the installed Julia packages and C++ dependencies:

```sh
julia --startup-file=no -O1 scripts/export_gsn_circular_m_mode.jl 0.5 8
julia --startup-file=no -O1 scripts/export_gsn_circular_m_mode.jl 0.0 8
julia --startup-file=no -O1 scripts/export_gsn_circular_m_mode.jl 0.5 6 Data/gsn_circular_a0.5_m2_lmax6.csv
cmake -S . -B build
cmake --build build --target bondi_kerr_m_mode_check -j4
build/bondi_kerr_m_mode_check Data/gsn_circular_a0.5_m2.csv Data/gsn_kerr_m2_decay_n13.csv 13
build/bondi_kerr_m_mode_check Data/gsn_circular_a0.5_m2.csv Data/gsn_kerr_m2_decay_n17.csv 17
build/bondi_kerr_m_mode_check Data/gsn_circular_a0.5_m2.csv Data/gsn_kerr_m2_decay_n21.csv 21
build/bondi_kerr_m_mode_check Data/gsn_circular_a0.0_m2.csv Data/gsn_schwarzschild_m2_decay.csv 13
build/bondi_kerr_m_mode_check Data/gsn_circular_a0.5_m2_lmax6.csv Data/gsn_kerr_m2_decay_lmax6.csv 13
ctest --test-dir build -R 'bondi_(kerr|schwarzschild_data|metric_falloff)' --output-on-failure
python3 scripts/plot_gsn_kerr_decay.py
```

The plotting script requires NumPy and Matplotlib. Outputs are
`plots/gsn_kerr/radial_decay.png` and `.pdf`. The Julia exporter currently fixes
the orbit radius and m pair as above; its command-line arguments set a, source
ell_max, and output path. The C++ executable requires odd Nz >=9 to include z=0.
