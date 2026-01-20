//
// Created by Peter Zimmerman on 25.10.25.
//

#ifndef GHZ_NUMERIC_GHPFIELDVECTORIZED_HPP
#define GHZ_NUMERIC_GHPFIELDVECTORIZED_HPP

// GHPFieldVectorized.hpp
// Modernized GHPField for vectorization

#pragma once

#include "SpinCoeffsNP.hpp"
#include "MathMacros.hpp"
#include "GHPScalars.hpp"

#include <vector>
#include <span>
#include <algorithm>
#include <execution>
#include <functional>
#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <concepts>
#include <type_traits>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace math;
/*!
 * \namespace ghp
 * \brief Geometrical-Held-Penrose (GHP) field container with vectorized storage
 *
 * GHPFieldVectorized stores a 2D field of GHP scalars in a flat 1D std::vector
 * for cache-friendly, contiguous memory layout. It supports element-wise
 * operations (addition, multiplication) and common transformations such as
 * conjugation and spin-boost transformations.
 *
 * \note OpenMP is used to parallelize loops over the 2D field for:
 *       - element-wise operations (operator+=, operator*=)
 *       - filling the field from indexed or nodal functions (fill_indexed, fill_nodes)
 *       - conjugation and spin-boost transformations
 *
 * \details
 * - `Nr_` and `Nz_` are the number of points in the radial and angular directions.
 * - `p_` and `q_` are the GHP weights of the field.
 * - Internally, the 2D array is flattened using `idx(r,z) = r*Nz_ + z`.
 * - All loops that touch every element of the field use `#pragma omp parallel for`
 *   with collapse(2) for 2D loops or a single loop for 1D.
 *
 * Usage example:
 * \code
 * ghp::GHPFieldVectorized field(Nr, Nz, Scalar(0.0,0,0), p, q);
 * field.fill_indexed([](int r, int z){ return Complex(r+z*I); });
 * field += other_field;
 * \endcode
 */
;
namespace ghp {

    using Real = teuk::Real; // replace with your preferred real type

    class GHPFieldVectorized {
    public:
        using Complex = teuk::Complex;
        using Scalar  = GHPScalar<Complex>;
        using container_type = std::vector<Scalar>;
        using value_type = Scalar;
        using size_type  = std::size_t;

    private:
        container_type values_;
        int p_ = 0;
        int q_ = 0;
        int Nz_ = 0;
        int Nr_ = 0;

        constexpr size_type idx(int r, int z) const noexcept {
            return static_cast<size_type>(r) * static_cast<size_type>(Nz_) + static_cast<size_type>(z);
        }

    public:
        GHPFieldVectorized() = default;

        explicit GHPFieldVectorized(int Nr, int Nz,
                                    Scalar init = Scalar(teuk::zeroC, 0, 0),
                                    int p = 0, int q = 0)
                : values_(static_cast<size_type>(Nr) * static_cast<size_type>(Nz)),
                  p_(p), q_(q),
                  Nz_(Nz), Nr_(Nr)
        {
            for (auto &s : values_) s = Scalar(init.value(), p_, q_);
        }

        int p() const noexcept { return p_; }
        int q() const noexcept { return q_; }
        int Nz() const noexcept { return Nz_; }
        int Nr() const noexcept { return Nr_; }

        Scalar& operator()(int r, int z) { return values_[idx(r, z)]; }
        const Scalar& operator()(int r, int z) const { return values_[idx(r, z)]; }

        void resize(int Nr, int Nz, Scalar init = Scalar(teuk::zeroC, 0, 0)) {
            Nr_ = Nr; Nz_ = Nz;
            values_.assign(static_cast<size_type>(Nr_) * static_cast<size_type>(Nz_),
                           Scalar(init.value(), p_, q_));
        }

/**
 * @brief Fill all (r,z) grid points by evaluating a user-supplied function.
 *
 * This routine loops over the full Nr_ × Nz_ grid and assigns each entry
 *     values_[idx(ir, iz)] using the value returned by the callable `f(ir, iz)`.
 *
 * The callable `f` must have the signature:
 *     teuk::Complex f(int ir, int iz)
 * or anything implicitly convertible to the value used by the Scalar constructor.
 *
 * Parallelism:
 *   The fill is executed with OpenMP over both r and z indices.
 *
 * Template parameters:
 *   F  - any callable satisfying std::invocable<int,int>
 *
 * Example usage:
 *
 * @code
 * HeldField psi(Nr, Nz, Scalar(0.0, p, q));
 *
 * psi.fill_indexed([&](int ir, int iz) {
 *     double r = r_nodes[ir];   // user-supplied nodes
 *     double z = z_nodes[iz];   // typically in [-1, 1]
 *
 *     // Example field value on the grid:
 *     return std::sin(r) * std::exp(z);
 * });
 * @endcode
 *
 * More complex usage (building from a tetrad or metric):
 *
 * @code
 * psi.fill_indexed([&](int ir, int iz) {
 *     X.r = r_nodes[ir];
 *     X.z = z_nodes[iz];
 *     tetrad.build_tetrad(X);
 *     return tetrad.get_Psi2(X);
 * });
 * @endcode
 *
 * @note The returned value is passed to:
 *         values_[idx(ir, iz)] = Scalar(f(ir, iz), p_, q_);
 *       so ensure the return type is compatible with the Scalar constructor.
 */
        template <std::invocable<int,int> F>
        void fill_indexed(F&& f) {
#pragma omp parallel for collapse(2)
            for (int ir = 0; ir < Nr_; ++ir)
                for (int iz = 0; iz < Nz_; ++iz)
                    values_[idx(ir,iz)] = Scalar(std::invoke(f, ir, iz), p_, q_);
        }
/**
 * @brief Fill the field values by evaluating a function on physical
 *        (r, z) node coordinates.
 *
 * This routine assigns every grid entry
 *     values_[idx(ir, iz)].value()
 * using the result of a user-supplied callable
 *
 *      f(r_nodes[ir], z_nodes[iz])
 *
 * The callable `f` must accept two `Real` arguments (r, z) and return
 * a value convertible to the underlying scalar type (the `.value()` part
 * of `Scalar`).
 *
 * Requirements:
 *   - `r_nodes.size() == Nr_`
 *   - `z_nodes.size() == Nz_`
 *   Otherwise an exception is thrown.
 *
 * Parallelism:
 *   OpenMP parallelization over both r and z indices.
 *
 * Template parameter:
 *   F — any callable satisfying `std::invocable<Real, Real>`
 *
 * Example usage:
 *
 * @code
 * SpectralField psi(Nr, Nz, Scalar(0.0, p, q));
 *
 * psi.fill_nodes(
 *     [&](Real r, Real z) {
 *         return std::sin(r) * std::cos(z);   // f(r,z)
 *     },
 *     r_nodes,    // span or vector of size Nr
 *     z_nodes     // span or vector of size Nz
 * );
 * @endcode
 *
 * Example using a tetrad / metric:
 *
 * @code
 * psi.fill_nodes(
 *     [&](Real r, Real z) {
 *         Coord X;
 *         X.r = r;
 *         X.z = z;
 *         tetrad.build_tetrad(X);
 *         return tetrad.get_Psi4(X);    // any scalar quantity
 *     },
 *     r_nodes,
 *     z_nodes
 * );
 * @endcode
 *
 * @note This function bypasses the full Scalar constructor and assigns only:
 *         values_[idx(r,z)].value() = f(...)
 *       The spin-weights (p_, q_) remain unchanged.
 */
        template <std::invocable<Real, Real> F>
        void fill_nodes(F&& f, std::span<const Real> r_nodes, std::span<const Real> z_nodes) {
            if (r_nodes.size() != static_cast<size_type>(Nr_) || z_nodes.size() != static_cast<size_type>(Nz_))
                throw std::invalid_argument("node vectors have wrong size");
#pragma omp parallel for collapse(2)
            for (int r = 0; r < Nr_; ++r)
                for (int z = 0; z < Nz_; ++z)
                    values_[idx(r,z)].value() = GHPScalar(std::invoke(f, r_nodes[r], z_nodes[z]), p_,q_);
        }

        GHPFieldVectorized& operator+=(const GHPFieldVectorized& other) {
            if (Nr_ != other.Nr_ || Nz_ != other.Nz_ || p_ != other.p_ || q_ != other.q_)
                throw std::runtime_error("Dimension or weight mismatch in operator+=");
#pragma omp parallel for
            for (size_type i = 0; i < values_.size(); ++i)
                values_[i] = values_[i] + other.values_[i];
            return *this;
        }

        friend GHPFieldVectorized operator+(const GHPFieldVectorized& a, const GHPFieldVectorized& b) {
            GHPFieldVectorized out = a;
            out += b;
            return out;
        }

        GHPFieldVectorized& operator*=(const GHPFieldVectorized& other) {
            if (Nr_ != other.Nr_ || Nz_ != other.Nz_)
                throw std::runtime_error("Dimension mismatch in operator*=");
#pragma omp parallel for  default(none) shared(other)
            for (size_type i = 0; i < values_.size(); ++i)
                values_[i] = values_[i] * other.values_[i];
            p_ += other.p_;
            q_ += other.q_;
            return *this;
        }

        friend GHPFieldVectorized operator*(const GHPFieldVectorized& a, const GHPFieldVectorized& b) {
            GHPFieldVectorized out(a);
            out *= b;
            return out;
        }

        [[nodiscard]] GHPFieldVectorized conj() const {
            GHPFieldVectorized out(Nr_, Nz_, Scalar(teuk::zeroC, 0, 0), q_, p_);
#pragma omp parallel for
            for (size_type i = 0; i < values_.size(); ++i)
                out.values_[i] = values_[i].conj();
            return out;
        }

        [[nodiscard]] GHPFieldVectorized transform(const Complex& lambda) const {
            GHPFieldVectorized out(Nr_, Nz_, Scalar(teuk::zeroC, 0, 0), p_, q_);
            const auto lambda_p  = math::PowInt(lambda, p_);
            const auto lambda_qc = math::PowInt(std::conj(lambda), q_);
#pragma omp parallel for
            for (size_type i = 0; i < values_.size(); ++i)
                out.values_[i] = lambda_p * lambda_qc * values_[i];
            return out;
        }

        [[nodiscard]] std::string str(int ir = -1, int jz = -1) const {
            std::ostringstream oss;
            oss << "GHPFieldVectorized(p,q)=(" << p_ << "," << q_ << "), size=(" << Nr_ << "," << Nz_ << ")\n";
            if (ir >= 0 && jz >= 0) {
                const auto &val = (*this)(ir, jz).value();
                oss << "values_[" << ir << "][" << jz << "] = " << std::setprecision(6) << val << '\n';
            }
            return oss.str();
        }
    };
// Example usage:
//
// ghp::GHPFieldVectorized A(100, 80, Scalar(teuk::zeroC,0,0), 1, -1);
// ghp::GHPFieldVectorized B(100, 80, Scalar(teuk::zeroC,0,0), 1, -1);
//
// // Fill with functions of index
// A.fill_indexed([](int r, int z){ return teuk::Complex(r + z, 0); });
// B.fill_indexed([](int r, int z){ return teuk::Complex(r - z, 0); });
//
// // Add
// auto C = A + B;
//
// // Multiply
// auto D = A * B;
//
// // Conjugate
// auto Ac = A.conj();
//
// // Spin-boost transform
// teuk::Complex lambda = {1.0, 0.5};
// auto AT = A.transform(lambda)
} // namespace ghp

namespace ghp {

    using Real = teuk::Real;
/*!
 * \class HeldFieldVectorized
 * \brief 1D vectorized GHP (Geometrical-Held-Penrose) field along the "z" direction
 *
 * HeldFieldVectorized is a lightweight 1D GHP field container that stores
 * scalars in a contiguous std::vector for cache-friendly access and
 * efficient parallel operations. It is conceptually a single-row
 * version of GHPFieldVectorized.
 *
 * \note OpenMP is used to accelerate operations over all elements:
 *       - fill_indexed()
 *       - fill_nodes()
 *       - operator+=, operator*=
 *       - conj()
 *
 * \details
 * - Nz_ stores the number of z-points.
 * - p_ and q_ store the GHP weights.
 * - Access elements via operator()(int z).
 * - Supports element-wise addition, multiplication, and conjugation.
 *
 * \example
 * ghp::HeldFieldVectorized h(10, ghp::GHPScalar(teuk::zeroC, 0, 0), 1, -1);
 *
 * // Fill using index
 * h.fill_indexed([](int z){ return teuk::Complex(z, z*0.1); });
 *
 * // Fill using node values
 * std::vector<ghp::Real> nodes = {0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9};
 * h.fill_nodes(nodes, [](ghp::Real z){ return teuk::Complex(z*z, z); });
 *
 * // Access element
 * auto val = h(3);
 * std::cout << h.str(3);
 */

    class HeldFieldVectorized {
    public:
        using Complex = teuk::Complex;
        using Scalar  = GHPScalar<Complex>;
        using container_type = std::vector<Scalar>;
        using value_type = Scalar;
        using size_type  = std::size_t;

    private:
        container_type values_;
        int p_ = 0;
        int q_ = 0;
        int Nz_ = 0;

    public:
        HeldFieldVectorized() = default;

        explicit HeldFieldVectorized(int Nz,
                                     Scalar init = Scalar(teuk::zeroC, 0, 0),
                                     int p = 0, int q = 0)
                : values_(static_cast<size_type>(Nz), Scalar(init.value(), p, q)),
                  p_(p), q_(q),
                  Nz_(Nz)
        {}

        [[nodiscard]] int p() const noexcept { return p_; }
        [[nodiscard]] int q() const noexcept { return q_; }
        [[nodiscard]] int Nz() const noexcept { return Nz_; }

        Scalar& operator()(int iz) { return values_[static_cast<size_type>(iz)]; }
        const Scalar& operator()(int iz) const { return values_[static_cast<size_type>(iz)]; }

        void resize(int Nz, Scalar init = Scalar(teuk::zeroC, 0, 0)) {
            Nz_ = Nz;
            values_.assign(static_cast<size_type>(Nz), Scalar(init.value(), p_, q_));
        }

        template <std::invocable<int> F>
        void fill_indexed(F&& f) {
#pragma omp parallel for
            for (int iz = 0; iz < Nz_; ++iz)
                values_[static_cast<size_type>(iz)] = Scalar(std::invoke(f, iz), p_, q_);
        }

        template <std::invocable<Real> F>
        void fill_nodes(std::span<const Real> z_nodes, F&& f) {
            if (z_nodes.size() != static_cast<size_type>(Nz_))
                throw std::invalid_argument("z_nodes vector has wrong size");
#pragma omp parallel for
            for (int iz = 0; iz < Nz_; ++iz)
                values_[static_cast<size_type>(iz)].value() = std::invoke(f, z_nodes[iz]);
        }

        HeldFieldVectorized& operator+=(const HeldFieldVectorized& other) {
            if (Nz_ != other.Nz_ || p_ != other.p_ || q_ != other.q_)
                throw std::runtime_error("Dimension or weight mismatch in operator+=");
#pragma omp parallel for
            for (size_type i = 0; i < values_.size(); ++i)
                values_[i] = values_[i] + other.values_[i];
            return *this;
        }

        friend HeldFieldVectorized operator+(const HeldFieldVectorized& a,
                                             const HeldFieldVectorized& b) {
            HeldFieldVectorized out = a;
            out += b;
            return out;
        }

        HeldFieldVectorized& operator*=(const HeldFieldVectorized& other) {
            if (Nz_ != other.Nz_)
                throw std::runtime_error("Dimension mismatch in operator*=");
#pragma omp parallel for
            for (size_type i = 0; i < values_.size(); ++i)
                values_[i] = values_[i] * other.values_[i];
            p_ += other.p_;
            q_ += other.q_;
            return *this;
        }

        friend HeldFieldVectorized operator*(const HeldFieldVectorized& a,
                                             const HeldFieldVectorized& b) {
            HeldFieldVectorized out(a);
            out *= b;
            return out;
        }

        [[nodiscard]] HeldFieldVectorized conj() const {
            HeldFieldVectorized out(Nz_, Scalar(teuk::zeroC, 0, 0), q_, p_);
#pragma omp parallel for
            for (size_type i = 0; i < values_.size(); ++i)
                out.values_[i] = values_[i].conj();
            return out;
        }

        std::string str(int iz = -1) const {
            std::ostringstream oss;
            oss << "HeldFieldVectorized(p,q)=(" << p_ << "," << q_ << "), size=" << Nz_ << "\n";
            if (iz >= 0) {
                const auto& val = values_[static_cast<size_type>(iz)].value();
                oss << "values_[" << iz << "] = " << std::setprecision(6) << val << '\n';
            }
            return oss.str();
        }

        const container_type& data() const { return values_; }
    };

} // namespace ghp
#endif //GHZ_NUMERIC_GHPFIELDVECTORIZED_HPP