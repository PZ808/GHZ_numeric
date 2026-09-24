# GeneralizedSasakiNakamura → Bondi reconstruction audit

Inspected 2026-09-22; mode-level numerical check added 2026-09-23.
This is an interface analysis with a Schwarzschild benchmark, not a complete converter.

## Numerical test, 2026-09-23

Ran GeneralizedSasakiNakamura 0.9.0 for M=mu=1, a=0, p=10, e=0,
x=1, (ell,m,n,k)=(2,+/-2,0,0), using `isem_trapezoidal`.
The installed package's complete `src/` tree compared identical to the
local checkout with `diff -qr` before the run.

For m=2:

    omega = 0.06324555320336758
    lambda = 4
    Zinf = -0.00112064079192542346 + 0.000305760838458181601 i

The independently evaluated m=-2 amplitude agrees with conjugate(Zinf)
to relative error 5.77e-16. Setting psi4^(1o)=-Zinf and applying the
direct equations below (A=B=6 for this Schwarzschild harmonic) gives

    f_direct    = 0.0017964781019034345 + 0.005849485943278971 i
    f_reference = 0.0017964735616899006 + 0.005849543375750201 i
    relative complex error = 9.414915918e-6.

The fbar coefficient is the same for this particular mode pair and agrees
to the same accuracy. No fitted normalization or phase was applied.
The equivalent psi0 coefficient inferred through the existing seed relation is

    inferred:  416.61441946399395 - 127.94947945199000 i
    CSV:       416.61850994235255 - 127.94915608709198 i.

This establishes agreement at approximately 1e-5 for one mode pair, not
machine precision or a generic-Kerr validation. The reference CSV is a
finite-radius fit; its printed precision is not an error bound. This test
alone does not identify the source of the remaining discrepancy or resolve
all signature versus rho-expansion interpretations of the reference.

The legacy `trapezoidal` calculation also completed successfully and agrees
with the ISEM complex amplitude to 5.54e-11 relative error for both signs
of m. It reproduces the approximately 9.41e-6 reference discrepancy.
This makes a solver-route-specific error unlikely; investigating the
reference fit/boundary accuracy is the next precision check. All assertions
in the reproduction script passed.

Reproduction: `scripts/check_gsn_schwarzschild_mode.jl` computes the pair,
compares both complex seeds, checks the original coupled equation and
also exercises the legacy `trapezoidal` route. Raw output is written to
`Data/gsn_schwarzschild_mode_check.csv`. This is a mode-level seed comparison;
it does not yet implement a new CSV importer into the C++ reconstruction.

## Correction: direct psi4 route (preferred for nonstatic modes)

The intermediate conversion to psi0 is optional. The user's supplied
relations, also appearing in Hollands–Toomani, *Metric reconstruction in
Kerr spacetime*, Appendix J and Appendix D, give a direct seed equation.
Source: https://inspirehep.net/files/636bc5854ae509d23b7da9f61024d163

Use the paper's mostly-minus conventions and let P=thorn-prime-Held,
A=eth-prime-Held^4 (spin +2 to -2), B=eth-Held^4 (spin -2 to +2).
The equations are

    A f - 3 Psi^o P fbar = P h_mbar_mbar^(1o)
    psi4^(1o) = (1/2) P^2 h_mbar_mbar^(1o).

For a nonzero-frequency mode exp(-i omega u + i m phi_*), P=-i omega
on these radially constant fields. Table B2 gives Psi^o=M. Hence

    A f + 3 i M omega fbar = (2 i / omega) psi4^(1o)
    B fbar + 3 i M omega f = (2 i / omega) conjugate(psi4)^(1o).

Here the conjugate field's coefficient at (m,omega) comes from the
conjugate of the original field at (-m,-omega), not conjugation at fixed
mode. For orbit labels the reflection is (m,n,k)->(-m,-n,-k); include
the chosen angular-basis conjugation phase when working with coefficients.
Eliminating f gives

    (A B + 9 M^2 omega^2) fbar
      = (2 i / omega) A conjugate(psi4)^(1o) + 6 M psi4^(1o).

Equivalently, eliminate fbar and use the existing spin +2 operator chain:

    (B A + 9 M^2 omega^2) f
      = (2 i / omega) B psi4^(1o) + 6 M conjugate(psi4)^(1o).

Thus the existing seed matrix can be reused with a different RHS, provided
the GHP weights and input conventions are matched. The modal reduction
has now been checked above; the generic angular-grid RHS is not implemented.

For the usual spin -2 Teukolsky master field T_-2=rho^(-4) psi4,
R_-2^out ~ Zinf r^3 exp(i omega r_*), and
rho^4=r^(-4)[1+O(1/r)]. Consequently the coefficient of 1/r in psi4
at fixed outgoing u is Zinf times the spin -2 spheroidal harmonic,
within the same Weyl/tetrad/sign and coordinate-phase conventions.

Crucially, the paper's expansion (50) uses powers of rho, not 1/r:

    psi4 = rho psi4^(1o) + O(rho^2)
         = -psi4^(1o)/r + O(r^(-2)).

The paper's psi4^(1o) is therefore minus the literal 1/r coefficient.
This odd-power sign is independent of a signature or Weyl-definition
conversion. Likewise distinguish a rho^5 coefficient from an r^(-5)
coefficient for psi0. Do not apply the WXF export's signature sign again
without independently establishing what its stored coefficient means.

Revised pipeline: Julia Zinf -> physical 1/r psi4 coefficient -> paper's
rho coefficient -> direct paired seed solve -> Bondi reconstruction.
Keep the psi0/TS route below as an independent validation path. The
static-sector and generic-Kerr metric limitations below still apply.

## Producer

`../tools/GeneralizedSasakiNakamura.jl/runmode.jl` calls
`Teukolsky_pointparticle_mode(-2,l,m,n,k,a,p,e,x)` and exports
`l,m,n,k,omega,lambda,Zinf_re,Zinf_im`. The public implementation is in
`src/GeneralizedSasakiNakamura.jl:1621`; `info.amplitude` is already a
Teukolsky amplitude, not a GSN amplitude requiring another conversion.
The spin +2 point-particle API returns the **horizon** amplitude, not the
spin +2 infinity amplitude required here.

The package uses M=1. The current script selects a=0.9, p=8, e=0.2,
x=cos(pi/4), ell<=8, |n|,|k|<=10: 33,957 mode requests. Its `scripts/`
directory is empty. `method="auto"` currently selects `isem_trapezoidal`
for this API. Homogeneous radial solutions alone do not supply the
point-particle normalization; the convolution integrals do.

`src/Homogeneous/ConversionFactors.jl:154` supplies
`TeukolskyStarobinsky_abs_Csq`, the radial TS modulus-squared quantity.
For s=+2 it shifts lambda by +4 to use the s=-2 convention. This real
quantity must not be confused with a complex TS conversion factor or
with the notebook's outgoing amplitude named Cplus.

## Consumer and notation

In `MathematicaNotebooks/EffectiveSource/lm_mode/BondiGauge_Adjusted_Generic.nb`:

- line 1440 sets `Weylspin=2`;
- around line 7551, `computeCplusMode` reads the Mathematica point-particle
  `Amplitudes[ScriptCapitalI]` for that spin;
- around line 9450, `Psi05HS[l,m,n,k]=-Cplus[Weylspin,l,m,n,k]`;
- around line 10490, `fHS=32 i omega^3 Psi05HS/StarobinskyScriptCapitalC`;
- `fbarHS` uses `(-1)^m Conjugate[fHS[L,-m,-n,-k,a]]`.

The nearby section heading mentions spin -2, but the executed definitions
use spin +2. Thus Zinf cannot simply be inserted into this Cplus cache.
The required interface is the physical TS conversion from spin -2 infinity
data to spin +2 infinity data, followed by the notebook's convention map.
Its sign, angular normalization, reflected-mode terms and phase conventions
remain to be validated; no scalar conversion factor is asserted here.
Radial TS identities between homogeneous basis functions alone do not fix
all of those physical-field conventions.

## Existing C++ entry point

`BondiHeldSeedSolver::solve(reduced_psi0,m,omega)` takes the pole-factorized
spin +2 angular field on the LGL nodes and solves

    [eth_H^4 ethbar_H^4 + 9 M^2 omega^2] f = 2 i omega^3 psi0^(5o).

The rightmost angular operator acts first. The Schwarzschild reference is

    D_l = (l-1) l (l+1) (l+2)
    f_lm = 32 i omega^3 psi_lm / [D_l^2 + 144 M^2 omega^2].

This is implemented by `SchwarzschildPsi0Modes::exact_reduced_seed`.
The package's TS modulus-squared reduces to the same denominator at a=0.

For Kerr, assemble each fixed (m,n,k) frequency separately:

    psi_reduced(z) = sum_l psi_lmnk * _2S_lm(a*omega;z) / P_m,2(z)
    P_m,s(z) = (1-z)^(|m+s|/2) (1+z)^(|m-s|/2).

Evaluate the reduced harmonics regularly at the poles; numerical division
of two vanishing quantities at an LGL endpoint is not acceptable. Sum over
ell at fixed frequency before the angular solve. Do not combine different
n,k frequencies into one fixed-m solve. Build the conjugate sector from
the reflected frequency mode, with the appropriate spin -2 basis.

`SchwarzschildPsi0Modes::load_csv` only accepts a complete (ell,m) table
with four columns and mostly-minus metadata. It cannot read the generic
Julia CSV. `BondiHeldMetricReconstruction` is explicitly a Schwarzschild
implementation, despite accepting Held operators; its current formulas
must not be assumed to implement the generic Kerr notebook. It also
rejects omega=0.

## Export issues to fix before production use

1. Resume checks only parse the first four columns. A truncated row with
   a valid key is incorrectly treated as complete. Validate all fields,
   finite values, duplicates, and orbit/convention metadata.
2. Orbit parameters appear to six decimal places in the filename. Use
   full metadata validation or a configuration hash to avoid collisions.
3. Record solver controls, package revision/path, units, particle-mass
   normalization, angular/tetrad/sign conventions, and orbit phase origins.
   Running the script without selecting the local Julia project can load
   the installed package rather than this checkout.
4. Keep successful, failed, symmetry-forbidden and static modes distinct.
   The infinity radiative solver skips |omega|<1e-12 and returns zero.
   This does not mean the static psi0 coefficient vanishes: the existing
   r0=10 table has psi_(2,0) approximately -342.53337749.
5. Preserve complex amplitudes and signed frequencies. Fluxes discard the
   phase needed for reconstruction. Do not enable low-flux zeroing without
   assessing errors in the converted field.

## Validation and implementation sequence

Start with a=0, p=10, e=0, x=1 and n=k=0. Compare TS-converted modes
(2,+/-2), (2,+/-1), and several ell=3 modes against
`tests/data/psi0_schwarzschild_r0_10_lmax20_mostly_minus.csv`. This exercises
phase, parity, overall sign and normalization; matching fluxes is insufficient.
The existing table was obtained by finite-radius fitting, so use a tolerance
consistent with that fit rather than its printed decimal precision.

Then feed converted nonstatic coefficients through the existing seed and
metric checks. Handle the static sector separately. Only after this anchor
is established extend the transport format to (ell,m,n,k), export/evaluate
spin +2 spheroidal functions, test nonzero-a angular identities, and port
the generic metric formulas from the notebook.

The existing built `bondi_schwarzschild_data_collocation_tests` executable
passed against the reference CSV during this audit: m=2 seed error
7.16e-10, residual 3.87e-8; m=2 and m=3 metric growing/constant terms
cancelled at approximately 2.15e-13 and 2.24e-12, with direct-sum radial
slopes close to -1. This verifies the existing reference-data path, not a
Julia-to-Bondi conversion, and was not a fresh rebuild.

A Julia (2,2,0,0) reference-mode calculation was started but interrupted
before returning an amplitude; no new Julia numerical result was obtained.
The default environment loaded the installed package under
`~/.julia/packages/GeneralizedSasakiNakamura/iL6U6`, rather than the checkout.
Its main module file compared identical to the checkout; the whole installed
source tree was not compared.
