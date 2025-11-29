//
// Created by Peter Zimmerman on 28.11.25.
//

#ifndef GHZ_NUMERIC_SPECTRALGHPFIELDVECTORIZED_HPP
#define GHZ_NUMERIC_SPECTRALGHPFIELDVECTORIZED_HPP


#pragma once

#include "GhzTypes.hpp"
#include "GHPScalars.hpp"
#include "HeldScalars.hpp"

#include <vector>
#include <array>
#include <span>
#include <functional>
#include <stdexcept>
#include <cassert>
#include <sstream>
#include <iomanip>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace ghp {

    using Real = teuk::Real;
    using Complex = teuk::Complex;

/**
 * @brief Generic 1D or 2D vectorized field container.
 *
 * FieldVectorized provides flat contiguous storage for scalars of type T,
 * with support for indexing, resizing, and parallel fill operations.
 *
 * Template Parameters:
 *   T   - the value type stored in the field (e.g., GHPScalar)
 *   Dim - dimensionality (1 or 2)
 */
    template<typename T, int Dim>
    class FieldVectorized {
    protected:
        std::vector<T> values_;
        std::array<int, Dim> dims_{};

        constexpr size_t linear_idx(const std::array<int, Dim>& idx) const noexcept {
            if constexpr (Dim == 1) return static_cast<size_t>(idx[0]);
            else if constexpr (Dim == 2) return static_cast<size_t>(idx[0]) * static_cast<size_t>(dims_[1]) + static_cast<size_t>(idx[1]);
            else static_assert(Dim <= 2, "FieldVectorized supports only 1D or 2D");
        }

    public:
        FieldVectorized() = default;

        FieldVectorized(const std::array<int, Dim>& dims, const T& init = T())
                : dims_(dims), values_(size_t(dims[0]) * (Dim==2 ? size_t(dims[1]) : 1), init) {}

        void resize(const std::array<int, Dim>& dims, const T& init = T()) {
            dims_ = dims;
            values_.assign(size_t(dims[0]) * (Dim==2 ? size_t(dims[1]) : 1), init);
        }

        const std::array<int, Dim>& dims() const noexcept { return dims_; }

        T& operator()(const std::array<int, Dim>& idx) { return values_[linear_idx(idx)]; }
        const T& operator()(const std::array<int, Dim>& idx) const { return values_[linear_idx(idx)]; }

        T* data_ptr_1d(int i) requires (Dim==2) { return &values_[i*dims_[1]]; }
        std::vector<T>& data() noexcept { return values_; }
        const std::vector<T>& data() const noexcept { return values_; }

        // fill with index-based callable
        template<std::invocable<int,int> F>
        void fill_indexed(F&& f) requires (Dim==2) {
#pragma omp parallel for collapse(2)
            for (int i = 0; i < dims_[0]; ++i)
                for (int j = 0; j < dims_[1]; ++j)
                    values_[linear_idx({i,j})] = std::invoke(f,i,j);
        }

        template<std::invocable<int> F>
        void fill_indexed(F&& f) requires (Dim==1) {
#pragma omp parallel for
            for (int i = 0; i < dims_[0]; ++i)
                values_[linear_idx({i})] = std::invoke(f,i);
        }

        // fill with node-based callable
        template<std::invocable<Real,Real> F>
        void fill_nodes(F&& f, std::span<const Real> x_nodes, std::span<const Real> y_nodes) requires (Dim==2) {
            if (x_nodes.size() != size_t(dims_[0]) || y_nodes.size() != size_t(dims_[1]))
                throw std::invalid_argument("Node vector size mismatch");
#pragma omp parallel for collapse(2)
            for (int i=0;i<dims_[0];++i)
                for (int j=0;j<dims_[1];++j)
                    values_[linear_idx({i,j})] = std::invoke(f,x_nodes[i],y_nodes[j]);
        }
    };

/**
 * @brief 1D vectorized GHP field along the "z" direction.
 *
 * HeldFieldVectorized stores GHP-scalars with spin weights (p,q) along
 * a single row. Supports element-wise arithmetic, conjugation, and
 * node-based filling. Flat contiguous storage enables cache-friendly access.
 */
    class HeldVectorized : public FieldVectorized<GHPScalar<Complex>,1> {
        int p_ = 0;
        int q_ = 0;
    public:
        using Base = FieldVectorized<GHPScalar<Complex>,1>;

        HeldVectorized() = default;
        HeldVectorized(int Nz,
                       GHPScalar<Complex> init=GHPScalar<Complex>(teuk::zeroC,0,0),
                               int p=0,int q=0)
                : Base({Nz}, init), p_(p), q_(q) {}

        int p() const noexcept { return p_; }
        int q() const noexcept { return q_; }

        GHPScalar<Complex>& operator()(int z) { return Base::operator()({z}); }
        const GHPScalar<Complex>& operator()(int z) const { return Base::operator()({z}); }

        // arithmetic operators
        HeldVectorized operator+(const HeldVectorized& other) const {
            if (p_ != other.p_ || q_ != other.q_ || dims_[0] != other.dims_[0])
                throw std::runtime_error("Dimension or spin mismatch");
            HeldVectorized out(dims_[0], GHPScalar<Complex>(teuk::zeroC,p_,q_), p_,q_);
#pragma omp parallel for
            for (int i=0;i<dims_[0];++i)
                out(i) = (*this)(i) + other(i);
            return out;
        }

    HeldVectorized operator*(const HeldVectorized& other) const {
            if (dims_[0] != other.dims_[0]) throw std::runtime_error("Dimension mismatch");
        HeldVectorized out(dims_[0], GHPScalar<Complex>(teuk::zeroC,0,0), p_+other.p_, q_+other.q_);
#pragma omp parallel for
            for (int i=0;i<dims_[0];++i)
                out(i) = (*this)(i) * other(i);
            return out;
        }

    HeldVectorized conj() const {
        HeldVectorized out(dims_[0], GHPScalar<Complex>(teuk::zeroC,0,0), q_,p_);
#pragma omp parallel for
            for (int i=0;i<dims_[0];++i)
                out(i) = (*this)(i).conj();
            return out;
        }
    };

/**
 * @brief 2D vectorized GHP field container.
 *
 * GHPFieldVectorized stores a 2D array of GHP-scalars with spin weights
 * (p,q) in a contiguous vector. Supports:
 *   - row-wise slicing (ZSlice) via std::span
 *   - element-wise arithmetic and conjugation
 *   - spin-boost transformations
 *   - OpenMP-parallelized fill operations
 */
    class GHPVectorized : public FieldVectorized<GHPScalar<Complex>,2> {
        int p_ = 0;
        int q_ = 0;
    public:
        using Base = FieldVectorized<GHPScalar<Complex>,2>;

        GHPVectorized() = default;
        GHPVectorized(int Nr, int Nz, GHPScalar<Complex> init=GHPScalar<Complex>(teuk::zeroC, 0, 0), int p=0, int q=0)
                : Base({Nr,Nz}, init), p_(p), q_(q) {}

        int p() const noexcept { return p_; }
        int q() const noexcept { return q_; }

        GHPScalar<Complex>& operator()(int r,int z) { return Base::operator()({r,z}); }
        const GHPScalar<Complex>& operator()(int r,int z) const { return Base::operator()({r,z}); }

        auto slice_z(int r) { return std::span<GHPScalar<Complex>>(Base::data_ptr_1d(r), dims_[1]); }

        // arithmetic operators
        GHPVectorized operator+(const GHPVectorized& other) const {
            if (p_ != other.p_ || q_ != other.q_ || dims_ != other.dims_) throw std::runtime_error("Dimension or spin mismatch");
            GHPVectorized out(dims_[0], dims_[1], GHPScalar<Complex>(teuk::zeroC,p_,q_), p_,q_);
#pragma omp parallel for collapse(2)
            for (int i=0;i<dims_[0];++i)
                for (int j=0;j<dims_[1];++j)
                    out(i,j) = (*this)(i,j) + other(i,j);
            return out;
        }

        GHPVectorized operator*(const GHPVectorized& other) const {
            if (dims_ != other.dims_) throw std::runtime_error("Dimension mismatch");
            GHPVectorized out(dims_[0], dims_[1], GHPScalar<Complex>(teuk::zeroC,0,0), p_+other.p_, q_+other.q_);
#pragma omp parallel for collapse(2)
            for (int i=0;i<dims_[0];++i)
                for (int j=0;j<dims_[1];++j)
                    out(i,j) = (*this)(i,j) * other(i,j);
            return out;
        }

        GHPVectorized conj() const {
            GHPVectorized out(dims_[0], dims_[1], GHPScalar<Complex>(teuk::zeroC,0,0), q_,p_);
#pragma omp parallel for collapse(2)
            for (int i=0;i<dims_[0];++i)
                for (int j=0;j<dims_[1];++j)
                    out(i,j) = (*this)(i,j).conj();
            return out;
        }
    };

/**
 * @brief 2D spectral field with mode metadata (w,m).
 *
 * SpectralFieldVectorized is a generic 2D field container for type T
 * (e.g., double or complex), augmented with mode metadata:
 *   - w: temporal/frequency mode
 *   - m: azimuthal mode
 *
 * Supports flat storage, slicing, resizing, and parallelized fill
 * operations.
 */
    template<typename T>
    class SpectralFieldVectorized : public FieldVectorized<T,2> {
    protected:
        int w_ = 0;
        int m_ = 0;
    public:
        using Base = FieldVectorized<T,2>;

        SpectralFieldVectorized() = default;
        SpectralFieldVectorized(int Nr,int Nz,int m,int w,const T& init=T())
                : Base({Nr,Nz}, init), w_(w), m_(m) {}

        int w() const noexcept { return w_; }
        int m() const noexcept { return m_; }

        T& operator()(int r,int z) { return Base::operator()({r,z}); }
        const T& operator()(int r,int z) const { return Base::operator()({r,z}); }

        auto slice_z(int r) { return std::span<T>(Base::data_ptr_1d(r), Base::dims()[1]); }
    };

/**
 * @brief 2D spectral GHP field with spin weights and mode metadata.
 *
 * SpectralGHPVectorized combines the features of GHPFieldVectorized
 * and SpectralFieldVectorized:
 *   - stores GHP-scalars with spin weights (p,q)
 *   - adds spectral mode metadata (w,m)
 *   - supports element-wise arithmetic and conjugation
 *   - provides row slicing via std::span
 *   - fully contiguous memory layout for cache efficiency
 */
    class SpectralGHPVectorized : public SpectralFieldVectorized<GHPScalar<Complex>> {
        int p_ = 0;
        int q_ = 0;
    public:
        using Base = SpectralFieldVectorized<GHPScalar<Complex>>;

        SpectralGHPVectorized() = default;
        SpectralGHPVectorized(int w,int m,int Nr, int Nz,
                              GHPScalar<Complex> init=GHPScalar<Complex>(
                                      teuk::zeroC,0,0), int p=0,int q=0)
                : Base(Nr,Nz,m,w,init), p_(p), q_(q) {}

        int p() const noexcept { return p_; }
        int q() const noexcept { return q_; }

        GHPScalar<Complex>& operator()(int r,int z) { return Base::operator()(r,z); }
        const GHPScalar<Complex>& operator()(int r,int z) const { return Base::operator()(r,z); }

        auto slice_z(int r) { return std::span<GHPScalar<Complex>>(Base::data_ptr_1d(r), Base::dims()[1]); }

        // arithmetic operators
        SpectralGHPVectorized operator+(const SpectralGHPVectorized& other) const {
            if (p_ != other.p_ || q_ != other.q_ || Base::dims() != other.Base::dims() || w_ != other.w() || m_ != other.m())
                throw std::runtime_error("Mismatch in spin, dims or modes");
            SpectralGHPVectorized out(Base::dims()[0], Base::dims()[1], m_, w_, GHPScalar<Complex>(teuk::zeroC,p_,q_), p_,q_);
#pragma omp parallel for collapse(2)
            for (int i=0;i<Base::dims()[0];++i)
                for (int j=0;j<Base::dims()[1];++j)
                    out(i,j) = (*this)(i,j) + other(i,j);
            return out;
        }

        SpectralGHPVectorized operator*(const SpectralGHPVectorized& other) const {
            if (Base::dims() != other.Base::dims()) throw std::runtime_error("Dimension mismatch");
            SpectralGHPVectorized out(Base::dims()[0], Base::dims()[1], m_, w_, GHPScalar<Complex>(teuk::zeroC,0,0), p_+other.p_,q_+other.q_);
#pragma omp parallel for collapse(2)
            for (int i=0;i<Base::dims()[0];++i)
                for (int j=0;j<Base::dims()[1];++j)
                    out(i,j) = (*this)(i,j) * other(i,j);
            return out;
        }

        SpectralGHPVectorized conj() const {
            SpectralGHPVectorized out(Base::dims()[0], Base::dims()[1], m_, w_, GHPScalar<Complex>(teuk::zeroC,0,0), q_,p_);
#pragma omp parallel for collapse(2)
            for (int i=0;i<Base::dims()[0];++i)
                for (int j=0;j<Base::dims()[1];++j)
                    out(i,j) = (*this)(i,j).conj();
            return out;
        }
    };

} // namespace ghp

#endif // GHZ_NUMERIC_SPECTRALGHPFIELDVECTORIZED_HPP