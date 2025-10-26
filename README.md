# GHZ_numeric

**Green Hollands Zimmerman (GHZ) Transport Solver**

A modular C++ framework for solving the **GHZ transport equations** 
for an effective source on Kerr spacetime.  
The code is designed for studies of **nonlinear gravitational perturbations**,
**self-force**, and **QNM mode coupling** 
using the nonlinear Teukolsky formalism and implements the 
**Geroch–Held–Penrose (GHP)**, and Held formalisms. 
---

## 🔭 Overview

`GHZ_numeric` provides analytic and numerical infrastructure for:
- Constructing the **Kerr metric** and its coordinate representations (Boyer–Lindquist, ingoing, and outgoing Kerr coordinates).
- Defining **null tetrads** (Kinnersley, Carter, Hartle–Hawking) and transforming between them.
- Computing **spin coefficients** in both NP and **GHP-covariant** form.
- Representing and manipulating **GHP scalars** (with spin- and boost-weights).
- Setting up the **transport equations** along null congruences (for the GHZ system).
- (In future modules) Importing the effective source in $m$-modes
- (In future modules) Integrating the GHZ shadowless **transport equations** numerically spectral and 
finite difference methods
- (In future modules) Reconstructing the effective metric from Teukolsky data and the corrector.
in a asymptotically flat shadowless radiation gauge 
- 
  The architecture is designed to be modular and extendable, allowing both **analytic** and **numerical** workflows.

---

## 🧱 Project Structure


Each class is independent and documented internally.  
The `main.cpp` file demonstrates how to:
- Construct the Kerr metric.
- Build a tetrad.
- Compute NP and GHP spin coefficients.
- Verify metric–tetrad consistency.

---

## ⚙️ Build Instructions

### Requirements
- **C++17** or newer (tested with Clang and GCC).
- **CMake ≥ 3.15**.
- (Optional) **Eigen3** for numerical linear algebra (future modules).

### Build

```bash
git clone https://github.com/<yourname>/GHZ_numeric.git
cd GHZ_numeric
mkdir build && cd build
cmake ..
make
