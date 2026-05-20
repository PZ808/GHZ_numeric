# Effective Teukolsky Source Projection

This README documents the current Mathematica workflow for computing spheroidal-harmonic projections of windowed effective Teukolsky sources from \(m\)-mode data.

The main output is a set of radial functions

\[
T^s_{\ell m}(X)
=
2\pi
\int_{-1}^{1} dy\,
\overline{{}_sS_{\ell m}(\arccos y,0;a\omega)}
\,
T^s_m(X,Y=-r_0 y),
\]

where

\[
y=\cos\theta,\qquad
Y=-r_0 y,\qquad
\omega=m\Omega_0.
\]

The current production run focuses on the \(s=+2\), \(m=2\) source with \(\ell=2,\ldots,15\), \(n_Y=128\) angular quadrature points, and 128 radial points per left/right patch.

---

## 1. Coordinate and notation conventions

The notebook uses three related but distinct variables:

```Mathematica
y      := Cos[theta]
Y      := -r0 y
X      := radial comoving coordinate (r-r0)/\mathcal{B}
```

Comoving parameters
\[
\rho_0 = r_0/M, \quad
\alpha = a/\Sqrt{Mr_0},\quad
\sqrt{\rho_0}\mathcal B = \sqrt{\alpha^2+\rho_0-2},\quad
\sqrt{\rho_0}\mathcal D = \sqrt{2\alpha+\rho_0-3}, \quad
z_c = 2\sqrt{M/r_0}\frac{\mathcal B}{\mathcal D \Omega_0}.
\]
In these parameters, we take a circular equatorial orbit
with 
\[
\Omega_0=\frac{1}{\sqrt{M\rho_0(\rho_0+\alpha)}}
\]


Important convention:

```Mathematica
yGrid  = Gauss-Legendre nodes for y = Cos[theta]
YGrid  = -r0 yGrid
```

The spheroidal harmonic is evaluated on `yGrid`, but the Teukolsky source function must be evaluated on `YGrid`.

Do **not** feed `yGrid` directly into `TTeuk...Compiled`.

The variable `z[]` used elsewhere in the local \(X,Y,z_c\) construction is not the same as \(y=\cos\theta\). It is a local/prolate-type variable and should not be used as the spheroidal projection variable.

Some references use \(y=\cos^2\theta\). That is **not** the convention in this notebook.

---

## 2. Source-building pipeline

The effective source is built in several stages.

### 2.1 Effective source components

The code constructs tetrad components of the effective source as \(m\)-mode functions of \((X,Y)\).

For spin \(-2\), the relevant components are

```Mathematica
{"nn", "nmb", "mbmb"}
```

For spin \(+2\), the relevant components are

```Mathematica
{"ll", "lm", "mm"}
```

### 2.2 Windowing

The local effective source must be windowed/localized before global angular projection. The unwindowed local extension is not guaranteed to be globally regular at the poles.

The current production source is windowed upstream, at the perturbation/effective-source construction level, not merely by multiplying the final Teukolsky source.

This distinction matters. Multiplying the final Teukolsky source by a window,

\[
T_s^{\rm Teuk}\to W T_s^{\rm Teuk},
\]

is only a numerical diagnostic. The consistent construction is to window before applying the differential Teukolsky source operator.

### 2.3 Teukolsky source expressions

After building the differential \(2\)-jets of the effective source components, the notebook constructs

```Mathematica
TTeukPlus2Expr[m]
TTeukMinus2Expr[m]
```

and takes a single test mode,

```Mathematica
mUse = Round[mModeval];

TTeukPlus2ExprFixed  = TTeukPlus2Expr[mUse];
TTeukMinus2ExprFixed = TTeukMinus2Expr[mUse];
```

The expressions are compiled for speed as

```Mathematica
TTeukPlus2Compiled
TTeukMinus2Compiled
```

with signatures

```Mathematica
TTeukPlus2Compiled[Xval, Yval]
TTeukMinus2Compiled[Xval, Yval]
```


---

## 3. Angular projection

The angular projection uses Gauss-Legendre quadrature on

\[
y=\cos\theta\in[-1,1].
\]

The angular projector is built by

```Mathematica
makeAngularProjector[s, ell, m, nY]
```

and returns an association containing

```Mathematica
<|
  "s" -> s,
  "ell" -> ell,
  "m" -> m,
  "omega" -> omega,
  "yGrid" -> yGrid,
  "YGrid" -> YGrid,
  "Weights" -> weights,
  "SpheroidalValues" -> sphVals,
  "Projector" -> projector
|>
```

where

```Mathematica
YGrid = -r0 yGrid
Projector = 2 Pi weights Conjugate[sphVals]
```

The projection at one radial point is

```Mathematica
projectAtXFast[Tfun, projector, Xval]
```

and over a radial grid:

```Mathematica
projectOnXGridFast[Tfun, projector, XGrid]
```

---

## 4. Radial patching

The radial coordinate is split into left and right patches around the particle radius:

```Mathematica
patchLeft["XGrid"]
patchRight["XGrid"]
```

For production runs with 128 radial points per side:

```Mathematica
nXLeftProd  = 128;
nXRightProd = 128;

radialGridsProd =
  makeRadialPatchGrids[
    nXLeftProd,
    nXRightProd,
    "EpsilonX" -> \[Epsilon]Tiny
  ];

XLeftGridProd  = radialGridsProd["XLeftGrid"];
XRightGridProd = radialGridsProd["XRightGrid"];
```

The production patch association is

```Mathematica
patchAssocProd =
 <|
  "Left"  -> patchLeftProd["XGrid"],
  "Right" -> patchRightProd["XGrid"]
 |>;
```

For clarity, avoid old/deprecated bare `XGrid` machinery. Use `patchLeft["XGrid"]`, `patchRight["XGrid"]`, or the production patch equivalents explicitly.

---

## 5. Parallel production loop for \(s=+2\)

The production loop parallelizes over coarse tasks,

\[
(\ell,\text{patch}),
\]

rather than over individual radial or angular points.

For \(\ell_{\max}=15\), \(m=2\), and two patches, there are

\[
14\times 2=28
\]

tasks.

Projectors are built on the master kernel and passed to subkernels. This avoids package-loading issues with `GaussianQuadratureWeights` on parallel kernels.

Run:

```Mathematica
ellMaxPlus = 15;
nYProd = 128;

plus2ProjectionResultsNX128 =
  runPlus2ProjectionsParallelWithPatches[
    ellMaxPlus,
    nYProd,
    patchAssocProd
  ];
```

The output is an association keyed by

```Mathematica
{+2, ell, "Left"}
{+2, ell, "Right"}
```

Example access:

```Mathematica
plus2ProjectionResultsNX128[{+2, 2, "Left"}]["Values"]
plus2ProjectionResultsNX128[{+2, 15, "Right"}]["Values"]
```

Each entry contains

```Mathematica
<|
  "s" -> +2,
  "ell" -> ell,
  "m" -> mUse,
  "Patch" -> patchName,
  "nY" -> nY,
  "XGrid" -> XGrid,
  "Values" -> vals,
  "TimeSeconds" -> timing,
  "SourceBuild" -> "windowed-at-perturbation-level"
|>
```

---

## 6. Diagnostics

### 6.1 Basic output check

```Mathematica
projectionEntryCheck[plus2ProjectionResultsNX128]
```

This checks:

- radial grid length,
- number of projected values,
- numerical validity,
- max/min amplitude.

### 6.2 Patch stitching near \(X=0\)

```Mathematica
Dataset[
 Table[
  patchEdgeCheck[resPlus, +2, ell],
  {ell, ellListPlus}
 ]
]
```

A healthy result has a small jump across the left/right split. In the current test, jumps were typically around

\[
10^{-5}\text{--}10^{-4},
\]

while projected values were around

\[
10^{-2}\text{--}10^{-1}.
\]

### 6.3 Plot individual modes

```Mathematica
plotProjectedMode[resPlus, +2, 2]
plotProjectedMode[resPlus, +2, 5]
plotProjectedMode[resPlus, +2, 10]
plotProjectedMode[resPlus, +2, 15]
```

The projected radial profiles should be smooth on each patch, localized near \(X=0\), and free of isolated outliers.

### 6.4 Plot all \(\ell\) modes

```Mathematica
ellListPlus = Range[Max[Abs[mUse], 2], 15];

absCurvesPlus =
  Table[
   getJoinedProjData[resPlus, +2, ell]["AbsData"],
   {ell, ellListPlus}
  ];

ListLinePlot[
 absCurvesPlus,
 PlotLegends -> Placed[ellListPlus, Right],
 PlotRange -> All,
 AxesLabel -> {"X", "|TTeukPlus2_{ell m}(X)|"}
]
```

### 6.5 Angular convergence

Before trusting production data, check representative modes against a higher angular reference:

```Mathematica
convPlus2Ell15Ref =
  projectionConvergenceToReference[
    TTeukPlus2Compiled,
    +2,
    15,
    mUse,
    patchLeftProd["XGrid"],
    {64, 96, 128, 160, 192}
  ];

convPlus2Ell15Ref["ErrorsToReference"]
```

The recommended production value is currently

```Mathematica
nY = 128
```

with spot checks against `nY = 192`.

---

## 7. Known pitfalls

### 7.1 `AssociationMap` on associations

Do not use

```Mathematica
AssociationMap[f, assoc]
```

when `assoc` is already an association. This can feed rules into functions unexpectedly.

Use

```Mathematica
Map[f, assoc]
```

or

```Mathematica
Association @ KeyValueMap[
  Function[{key, val}, key -> f[val]],
  assoc
]
```

### 7.2 `GaussianQuadratureWeights` on subkernels

Avoid building projectors on parallel subkernels. Build all projectors on the master kernel and include the completed projector association in each task.

This avoids errors like

```Mathematica
GaussianQuadratureWeights[...][[All,1]]
is longer than depth of object
```

which occur when the package is not loaded on a subkernel.

### 7.3 `Compile::part` from unresolved symbolic parts

Warnings of the form

```Mathematica
EFEmhS1MMode[[5]] cannot be compiled
```

mean that symbolic `Part` expressions survived into a compiled expression. Force component expressions to evaluate before compiling.

Diagnostic:

```Mathematica
Cases[
 Hold[teffJet2["nn", mUse]],
 HoldPattern[EFEmhS1MMode[[__]]],
 Infinity
] // DeleteDuplicates
```

This should return `{}` before compilation.

### 7.4 Endpoint behavior

The unwindowed local source extension produced bad polar behavior under global spheroidal projection. The projection only became well-behaved after windowing/localizing the source upstream.

The endpoint issue is not a quadrature bug: it reflects using a local extension as though it were globally regular on the full sphere.

---

## 8. Saving results


Recommended metadata to include in exported structures:

```Mathematica
<|
 "s" -> +2,
 "m" -> mUse,
 "ellRange" -> {2, 15},
 "nY" -> 128,
 "nXPerSide" -> 128,
 "SourceBuild" -> "windowed-at-perturbation-level",
 "CoordinateConvention" -> "y = Cos[theta], Y = -r0 y",
 "ProjectionFormula" -> "2 Pi int_{-1}^{1} dy Conjugate[S_s,l,m(y)] T_s,m(X,-r0 y)"
|>
```

---

## 9. Current status

The \(s=+2\), \(m=2\), \(\ell=2,\ldots,15\) projection with

```Mathematica
nY = 128
nX = 128 per side
```

looks qualitatively healthy:

- smooth radial profiles,
- localized near \(X=0\),
- no polar endpoint blow-up after windowing,
- small left/right patch jumps,
- suitable for a pilot handoff with convergence diagnostics.
                                       
