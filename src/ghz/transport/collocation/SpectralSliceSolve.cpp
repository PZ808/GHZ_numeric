//
// Created by Peter Zimmerman on 19.03.26.
//

#include "ghz/transport/collocation/SpectralSliceSolve.hpp"

#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <algorithm>
#include <sstream>
#include <stdexcept>

namespace ghz::collocation {

    namespace ublas = boost::numeric::ublas;

    namespace {

        using CMatrix = ublas::matrix<Complex>;
        using PMatrix = ublas::permutation_matrix<std::size_t>;

        // ---------------------------------------------------------------------
        // Solve A x = b for complex dense square A.
        // ---------------------------------------------------------------------
        std::vector<Complex> solve_complex_linear_system(
                const CMatrix& A,
                const std::vector<Complex>& b)
        {
            const std::size_t N = b.size();
            if (A.size1() != N || A.size2() != N) {
                throw std::runtime_error(
                        "solve_complex_linear_system: matrix/vector size mismatch.");
            }

            CMatrix LU(A);
            PMatrix pm(N);

            const int status = ublas::lu_factorize(LU, pm);
            if (status != 0) {
                throw std::runtime_error(
                        "solve_complex_linear_system: LU factorization failed.");
            }

            ublas::vector<Complex> x(N);
            for (std::size_t i = 0; i < N; ++i) {
                x(i) = b[i];
            }

            ublas::lu_substitute(LU, pm, x);

            std::vector<Complex> out(N);
            for (std::size_t i = 0; i < N; ++i) {
                out[i] = x(i);
            }
            return out;
        }

        // ---------------------------------------------------------------------
        // Consistency checks
        // ---------------------------------------------------------------------
        void check_scalar_equation_sizes(
                const ScalarEquationSlice& eq,
                std::size_t N,
                const char* where)
        {
            if (eq.a2.size() != N || eq.a1.size() != N ||
                eq.a0.size() != N || eq.rhs.size() != N) {
                std::ostringstream os;
                os << where << ": equation vector size mismatch with domain size "
                   << N;
                throw std::runtime_error(os.str());
            }
        }

        void check_domain_sizes(std::size_t NL, std::size_t NR)
        {
            // We overwrite:
            //  - left row 0  (outer BC 1)
            //  - left row 1  (outer BC 2)
            //  - left row NL-1 (interface value row)
            //  - right row 0 (interface derivative row)
            //
            // so we need at least NL >= 3, NR >= 1.
            // In practice NR >= 2 is a safer minimal collocation solve.
            if (NL < 3) {
                throw std::runtime_error(
                        "solve_scalar_two_domain_slice: left domain must have at least 3 nodes.");
            }
            if (NR < 2) {
                throw std::runtime_error(
                        "solve_scalar_two_domain_slice: right domain must have at least 2 nodes.");
            }
        }

        // ---------------------------------------------------------------------
        // Fill block-diagonal collocation matrix from the two scalar equations.
        // ---------------------------------------------------------------------
        void assemble_block_operator(
                CMatrix& B,
                std::vector<Complex>& rhs,
                const ghz::numeric::PhysicalChebRadialOps& rops_left,
                const ghz::numeric::PhysicalChebRadialOps& rops_right,
                const TwoDomainEquationSlice& eq)
        {
            const auto& DrL  = rops_left.Dr_matrix();
            const auto& DrrL = rops_left.Drr_matrix();
            const auto& DrR  = rops_right.Dr_matrix();
            const auto& DrrR = rops_right.Drr_matrix();

            const std::size_t NL = rops_left.r().size();
            const std::size_t NR = rops_right.r().size();

            // Left block rows: 0 .. NL-1
            for (std::size_t i = 0; i < NL; ++i) {
                rhs[i] = eq.left.rhs[i];

                for (std::size_t j = 0; j < NL; ++j) {
                    const Complex Iij = (i == j) ? Complex(1.0, 0.0) : Complex(0.0, 0.0);
                    B(i, j) = eq.left.a2[i] * Complex(DrrL(i, j), 0.0)
                              + eq.left.a1[i] * Complex(DrL(i, j), 0.0)
                              + eq.left.a0[i] * Iij;
                }
            }

            // Right block rows: NL .. NL+NR-1
            for (std::size_t i = 0; i < NR; ++i) {
                rhs[NL + i] = eq.right.rhs[i];

                for (std::size_t j = 0; j < NR; ++j) {
                    const Complex Iij = (i == j) ? Complex(1.0, 0.0) : Complex(0.0, 0.0);
                    B(NL + i, NL + j) = eq.right.a2[i] * Complex(DrrR(i, j), 0.0)
                                        + eq.right.a1[i] * Complex(DrR(i, j), 0.0)
                                        + eq.right.a0[i] * Iij;
                }
            }
        }

        // ---------------------------------------------------------------------
        // Overwrite one row with an outer boundary condition.
        // ---------------------------------------------------------------------
        void overwrite_outer_bc_row(
                CMatrix& B,
                std::vector<Complex>& rhs,
                const ghz::numeric::PhysicalChebRadialOps& rops_left,
                const ghz::numeric::PhysicalChebRadialOps& rops_right,
                const BoundaryCondition& bc,
                std::size_t row)
        {
            const auto& DrL = rops_left.Dr_matrix();
            const auto& DrR = rops_right.Dr_matrix();

            const std::size_t NL = rops_left.r().size();
            const std::size_t NR = rops_right.r().size();
            const std::size_t NT = NL + NR;

            for (std::size_t j = 0; j < NT; ++j) {
                B(row, j) = Complex(0.0, 0.0);
            }

            if (bc.side == BCSide::Left) {
                const std::size_t idxL = 0;

                if (bc.kind == BCKind::Value) {
                    B(row, idxL) = Complex(1.0, 0.0);
                } else {
                    for (std::size_t j = 0; j < NL; ++j) {
                        B(row, j) = Complex(DrL(idxL, j), 0.0);
                    }
                }
            } else {
                const std::size_t idxR = NR - 1;

                if (bc.kind == BCKind::Value) {
                    B(row, NL + idxR) = Complex(1.0, 0.0);
                } else {
                    for (std::size_t j = 0; j < NR; ++j) {
                        B(row, NL + j) = Complex(DrR(idxR, j), 0.0);
                    }
                }
            }

            rhs[row] = bc.value;
        }

        // ---------------------------------------------------------------------
        // Overwrite interface rows.
        //
        // We use:
        //   row (NL-1) : u_R(r_p)  - u_L(r_p)  = value_jump
        //   row (NL)   : u'_R(r_p) - u'_L(r_p) = derivative_jump
        //
        // Here:
        //   left interface node  = last node of left domain
        //   right interface node = first node of right domain
        // ---------------------------------------------------------------------
        void overwrite_interface_rows(
                CMatrix& B,
                std::vector<Complex>& rhs,
                const ghz::numeric::PhysicalChebRadialOps& rops_left,
                const ghz::numeric::PhysicalChebRadialOps& rops_right,
                const InterfaceCondition& iface)
        {
            const auto& DrL = rops_left.Dr_matrix();
            const auto& DrR = rops_right.Dr_matrix();

            const std::size_t NL = rops_left.r().size();
            const std::size_t NR = rops_right.r().size();
            const std::size_t NT = NL + NR;

            const std::size_t row_value = NL - 1;
            const std::size_t row_deriv = NL;

            const std::size_t iL = NL - 1; // last node on left
            const std::size_t iR = 0;      // first node on right

            // Value jump row
            for (std::size_t j = 0; j < NT; ++j) {
                B(row_value, j) = Complex(0.0, 0.0);
            }
            B(row_value, iL)     = Complex(-1.0, 0.0);
            B(row_value, NL + iR)= Complex( 1.0, 0.0);
            rhs[row_value] = iface.value_jump;

            // Derivative jump row
            for (std::size_t j = 0; j < NT; ++j) {
                B(row_deriv, j) = Complex(0.0, 0.0);
            }

            for (std::size_t j = 0; j < NL; ++j) {
                B(row_deriv, j) = Complex(-DrL(iL, j), 0.0);
            }
            for (std::size_t j = 0; j < NR; ++j) {
                B(row_deriv, NL + j) = Complex(DrR(iR, j), 0.0);
            }

            rhs[row_deriv] = iface.derivative_jump;
        }

    } // namespace

    // -------------------------------------------------------------------------
    // Solve one z-slice
    // -------------------------------------------------------------------------
    TwoDomainSolutionSlice solve_scalar_two_domain_slice(
            const ghz::numeric::PhysicalChebRadialOps& rops_left,
            const ghz::numeric::PhysicalChebRadialOps& rops_right,
            const TwoDomainEquationSlice& eq,
            const BoundaryCondition& bc_row0,
            const BoundaryCondition& bc_row1,
            const InterfaceCondition& iface)
    {
        const std::size_t NL = rops_left.r().size();
        const std::size_t NR = rops_right.r().size();
        const std::size_t NT = NL + NR;

        check_domain_sizes(NL, NR);
        check_scalar_equation_sizes(eq.left,  NL, "solve_scalar_two_domain_slice(left)");
        check_scalar_equation_sizes(eq.right, NR, "solve_scalar_two_domain_slice(right)");

        CMatrix B(NT, NT, Complex(0.0, 0.0));
        std::vector<Complex> rhs(NT, Complex(0.0, 0.0));

        assemble_block_operator(B, rhs, rops_left, rops_right, eq);

        // Replace rows:
        //   0      -> outer BC row 0
        //   1      -> outer BC row 1
        //   NL-1   -> interface value row
        //   NL     -> interface derivative row
        overwrite_outer_bc_row(B, rhs, rops_left, rops_right, bc_row0, 0);
        overwrite_outer_bc_row(B, rhs, rops_left, rops_right, bc_row1, 1);
        overwrite_interface_rows(B, rhs, rops_left, rops_right, iface);

        const auto u = solve_complex_linear_system(B, rhs);

        TwoDomainSolutionSlice sol;
        sol.left.resize(NL);
        sol.right.resize(NR);

        for (std::size_t i = 0; i < NL; ++i) {
            sol.left[i] = u[i];
        }
        for (std::size_t i = 0; i < NR; ++i) {
            sol.right[i] = u[NL + i];
        }

        return sol;
    }

    // -------------------------------------------------------------------------
    // Solve all z-slices
    // -------------------------------------------------------------------------
    TwoDomainSolutionGrid solve_scalar_two_domain_grid(
            const ghz::numeric::PhysicalChebRadialOps& rops_left,
            const ghz::numeric::PhysicalChebRadialOps& rops_right,
            size_t Nz,
            const Modes& modes,
            ghp::GHPType out_type,
            const TwoDomainSliceEquationBuilder& eq_builder,
            const SliceBCBuilder& bc_row0_builder,
            const SliceBCBuilder& bc_row1_builder,
            const InterfaceConditionBuilder& iface_builder)
    {
        if (!eq_builder) {
            throw std::runtime_error(
                    "solve_scalar_two_domain_grid: eq_builder not set.");
        }
        if (!bc_row0_builder || !bc_row1_builder) {
            throw std::runtime_error(
                    "solve_scalar_two_domain_grid: boundary condition builders not set.");
        }
        if (!iface_builder) {
            throw std::runtime_error(
                    "solve_scalar_two_domain_grid: interface condition builder not set.");
        }

        const std::size_t NL = rops_left.r().size();
        const std::size_t NR = rops_right.r().size();

        TwoDomainSolutionGrid out{
                SpectralField(
                        NL, Nz, modes,
                        ghp::GHPScalar<Complex>(teuk::zeroC, out_type.p, out_type.q),
                        out_type.p, out_type.q
                ),
                SpectralField(
                        NR, Nz, modes,
                        ghp::GHPScalar<Complex>(teuk::zeroC, out_type.p, out_type.q),
                        out_type.p, out_type.q
                )
        };

#pragma omp parallel for
        for (size_t iz = 0; iz < Nz; ++iz) {
            const auto eq    = eq_builder(iz);
            const auto bc0   = bc_row0_builder(iz);
            const auto bc1   = bc_row1_builder(iz);
            const auto iface = iface_builder(iz);

            const auto sol = solve_scalar_two_domain_slice(
                    rops_left, rops_right, eq, bc0, bc1, iface);

            for (std::size_t ir = 0; ir < NL; ++ir) {
                out.left.set_index(
                        ir, iz,
                        ghp::GHPScalar<Complex>(sol.left[ir], out_type.p, out_type.q)
                );
            }
            for (std::size_t ir = 0; ir < NR; ++ir) {
                out.right.set_index(
                        ir, iz,
                        ghp::GHPScalar<Complex>(sol.right[ir], out_type.p, out_type.q)
                );
            }
        }

        return out;
    }

    // -------------------------------------------------------------------------
    // Residual norm helper
    // -------------------------------------------------------------------------
    Real max_residual_inf_norm(
            const ublas::matrix<Complex>& B,
            const std::vector<Complex>& u,
            const std::vector<Complex>& rhs)
    {
        const std::size_t N = u.size();
        if (B.size1() != N || B.size2() != N || rhs.size() != N) {
            throw std::runtime_error(
                    "max_residual_inf_norm: size mismatch.");
        }

        Real max_res = Real(0);
        for (std::size_t i = 0; i < N; ++i) {
            Complex sum = Complex(0.0, 0.0);
            for (std::size_t j = 0; j < N; ++j) {
                sum += B(i, j) * u[j];
            }
            max_res = std::max(max_res, std::abs(sum - rhs[i]));
        }
        return max_res;
    }

} // namespace ghz::collocation