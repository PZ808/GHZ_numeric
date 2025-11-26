# GHZ_numeric

**Green Hollands Zimmerman (GHZ) Transport Solver**

A modular C++ framework for solving the **GHZ transport equations** 
for a generic orbit effective source in Kerr spacetime.  

---

## 🔭 Overview

`GHZ_numeric` provides numerical infrastructure for:
- Constructing the **Kerr metric** and its coordinate representations (Boyer–Lindquist, ingoing, and outgoing Kerr coordinates).
- Defining **null tetrads** (Kinnersley, Carter, Hartle–Hawking) and transforming between them.
- Computing **spin coefficients** in both NP and **GHP-covariant** form.
- Representing and manipulating **GHP scalars** (with spin- and boost-weights).
- Setting up the **transport equations** along null congruences (for the GHZ system).
- (In future modules) Importing the effective source in $m$-modes
- (In future modules) Integrating the GHZ **transport equations** numerically using spectral and 
finite difference methods
- (In future modules) Reconstructing the effective metric from Teukolsky data and the corrector.


  The architecture is designed to be modular and extendable, allowing both **analytic** and **numerical** workflows.

---

## 🧱 Project Structure

```mermaid
flowchart TD
    subgraph Geometry
        A[Metric] --> B[KerrMetric] 
        C[CoordinateSystem] --> D[KerrCharts]
        C --> E[KinnersleyTetrad]
    end
    subgraph Scalars
        E --> F[GHPScalar]
        F --> G[GHPField]
        G --> H[SpectralField] --> I[GHPSpectralField]
        J[HeldField]  --> K[SpectralDiffer]
    end
```  

Each class is independent and documented internally.  
The `main.cpp` file demonstrates how to:
- Construct the Kerr metric.
- Build a tetrad.
- Compute NP and GHP spin coefficients.
- Verify metric–tetrad consistency.

---


# Core Components

Below are the main classes and their purposes.

---

## 1. `KerrMetric`
Coordinate-independent Kerr metric object.

- Stores parameters \(M, a, \)
- Provides metric functions and invariants
- Backend used by all coordinate-specific metrics

---

## 2. `CoordinateSystem`
Constructs quantities in and provides functions to transforms between:

- Boyer–Lindquist coordinates
- Ingoing Kerr coordinates
- Outgoing Kerr coordinates
- Outgoing conformally compactified coordinates

Provides:

- Metric construction in each coordinate chart
- Transformation Jacobians
- Useful geometric quantities for tetrads and scalars

---

## 3. `KinnersleyTetrad<Coordtype T>?`
Builds the Kinnersley tetrad in all supported coordinates.

Computes:

- Newman–Penrose spin coefficients
- GHP spin/boost weights
- Held coefficients
- Basis vectors for acting with GHP operators

---

## 4. `GHPScalar<Complex T>`
Complex scalar with GHP weights \((p,q)\).

- Operator overloads: `+, -, *, /`
- Correct GHP transformation behavior
- Type-safe representation of weighted scalars

This is the basic algebraic object used everywhere.

---

## 5. `GHPField`a
Represents a field of 
Geroch-Held-Penrose (GHP) scalars with spin-boost weights (p,q)
on an \(N_r \times N_z\)   grid of \((r,z)\) values.

- Stored as a 2D grid `[r][z]`
- Arithmetic operator overloads which handle GHP weights 
- accessors and views
- lambda functions which fill the values given a function of r,z
- GHP covariant tranformations inc as conjugation and spin-boosts
- Used for background quantities such as GHP spin coefficients

---

## 6. `SpectralField<T>`
Generic container for spectral fields.

Important Definitions:

- `ZRow` = spectral points in \(z=\cos\theta\)
- `RZBlock` = 2D array `[r][z]`
- `ValueType = T` defines data type of the field (tuek::Complex for multiprecision)
- `w` = \(\omega\)
- `m` = \(m\)

Features:

- Supports arbitrary scalar type (e.g. GHPScalar, Real, Complex)
- Stores \((m, \omega)\) mode labels
- Provides access to slices, fill operations, and arithmetic

---

## 7. `GHPSpectralField`
A `SpectralField<GHPScalar>` with GHP bookkeeping.

- Tracks \((p,q)\) weights
- Tracks \((m,\omega)\) indices
- Supports GHP operator application in spectral form

---

## 8. `HeldField`
Held-scalar version of GHPField.

- Used for outgoing-ray ODEs and puncture subtraction
- Inherits GHP structure but with Held-weight rules

---

## 9. `SpectralDiffer`
Legendre–Gauss–Lobatto collocation, barycentric interpolation and differentiation.

Provides:

- LGL node construction in the z-direction
-  d/dz via Legendre differentiation matrices
- Edth and edth′ operators
- Barycentric interpolation

Used to operate on `ZSlice` objects of a `SpectralField`.

---

## 10. `KerrBoundOrbit`
Action–angle parametrization of bound orbits in Kerr spacetime.

- Computes mino and BL frequencies via elliptic integrals
- Decomposes motion into secular and oscillatory pieces 
- Builds Fourier decomposition of phase variables 
- Keplerian parametrization
- Supplies data needed for constructing puncture sources

---

## ⚙️ Build Instructions

### Requirements
- **C++17** or newer (tested with Clang and GCC).
- **CMake ≥ 3.15**.
-  **Boost** for multiprecision, elliptic integral, and linear algebra

### Build

```bash
git clone https://github.com/<yourname>/GHZ_numeric.git
cd GHZ_numeric
mkdir build && cd build
cmake ..
make

