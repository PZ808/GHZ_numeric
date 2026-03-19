//
// Created by Peter Zimmerman on 28.11.25.
//
#ifndef GHZ_NUMERIC_SPECTRALGHPFIELDVECTORIZED_HPP
#define GHZ_NUMERIC_SPECTRALGHPFIELDVECTORIZED_HPP


#pragma once

#include "ghz/core/GhzTypes.hpp"
#include "ghz/ghp/GHPScalars.hpp"
#include "ghz/ghp/HeldScalars.hpp"
#include "ghz/orbit/SourceTorusFrequencies.hpp"

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

namespace spectral {

    using Real = teuk::Real;
    using Complex = teuk::Complex;
    using namespace teuk::literals;
    using ghp::GHPScalar;

    // global time-stepping parameters for evolution tests
    const int Nt = 10;
    const Real dt = 0.1_r;
    const Real total_time = Nt*dt;
    const Real f_s = 1.0_r/total_time; // sampling frequency
    const Real f_N = 0.5_r/total_time; // Nyquist frequency

/**
 * @brief Generic 1D or 2D vectorized field container.
 * @class FieldVectorized
 *
 * FieldVectorized provides flat contiguous storage for scalars of type T,
 * with support for indexing, resizing, and parallel fill operations.
 *
 * @param
 * dims - array specifying the size in each dimension
 * @param
 * init - initial value to fill the field with
 *
 * @tparam
 *   T   - the value type stored in the field (e.g., GHPScalar)
 *  @tparam
 *   Dim - dimensionality (1 or 2)
 */
    template<typename T, size_t Dim>
    class FieldVectorized {

    protected:
        std::vector<T> values_;
        std::array<size_t, Dim> dims_{};

/**
 *
 * @brief Convert a multi-dimensional index to a linear index for vector storage.
 * @function linear_idx
 *
 * This helper maps a 1D or 2D integer index into the corresponding
 * flat-array offset used by FieldVectorized's internal storage.
 *
 * say we have a 2D field with Nz=4, Nr=4: dims_ = {4,4} (Nr rows, Nz cols)
 * Row 0:  0   1   2   3
 * Row 1:  4   5   6   7
 * Row 2:  8   9  10  11
 * Row 3: 12  13  15  16
 * this gets stores in row-major order as a flat array:
 * [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
 * so the linear index mapping is for Dim=2: is linear_idx({i,j}) = i * dims_[1] + j
 * Here are the mappings:
 *  idx (i,j)   linear_idx
 *  (0,0)  ->    0
 *  (0,1)  ->    1
 *  (0,2)  ->    2
 *  (0,3)  ->    3
 *  (1,0)  ->    4
 *  (1,1)  ->    5
 *  (1,2)  ->    6
 *  (1,3)  ->    7
 *  (2,0)  ->    8
 *      ...
 * (i,j)  ->   i * dims_[1] + j
 *
 *
 * For Dim = 1:
 *      linear_idx({i}) = i
 *
 * For Dim = 2:
 *      linear_idx({i, j}) = i * dims_[1] + j
 *
 * The function is evaluated at compile time whenever possible
 * (constexpr), and will trigger a static assertion if Dim > 2.
 * dims_[1] is the size of the second dimension (e.g., Nz for a 2D field).
 * @param idx  Multi-dimensional index (size Dim).
 * @return     Linearized index into the underlying 1D storage.
 */
        [[nodiscard]] constexpr size_t linear_idx(const std::array<size_t, Dim> &idx) const noexcept {
            if constexpr (Dim == 1) return static_cast<size_t>(idx[0]);
            else if constexpr (Dim == 2)
                return static_cast<size_t>(idx[0]) * static_cast<size_t>(dims_[1])
                       + static_cast<size_t>(idx[1]);
            else static_assert(Dim <= 2, "FieldVectorized supports only 1D or 2D");
        }

    public:
        FieldVectorized() = default;
        FieldVectorized(const std::array<size_t, Dim> &dims, const T &init = T())
                : dims_(dims), values_(size_t(dims[0]) * (Dim==2 ? size_t(dims[1]) : 1), init) {}

        ~FieldVectorized() = default; // Default destructor

        // resize the field to new dimensions and initialize
        void resize(const std::array<size_t, Dim> &dims, const T &init = T()) {
            dims_ = dims;
            values_.assign( size_t(dims[0]) * (Dim==2 ? size_t(dims[1]):1), init );
        }

        // metadata: dims
        [[nodiscard]] const std::array<size_t, Dim> &dims() const noexcept { return dims_; }  // return dimensions


        //
        // data access
        //

        // index access operators – must pass array of indices {..}
        T &operator()(const std::array<size_t, Dim> &idx) {
            return values_[linear_idx(idx)];
        }  // must pass array of indices {..}
        const T &operator()(const std::array<size_t, Dim> &idx) const { return values_[linear_idx(idx)]; }

        // pointer to contiguous row (fixed r)
        // Get a pointer to row i in a 2D array stored as a 1D contiguous block (row major = row i starts at i * dims_[1])
        T* data_ptr_1d(size_t i) requires (Dim == 2) { return &values_[i * dims_[1]]; }  // pointer to row i
        T* data_ptr_1d(int i) requires  (Dim == 2) { return &values_[i * dims_[1]]; }    // pointer to row i
        [[nodiscard]] const T* data_ptr_1d(size_t i) const requires (Dim == 2) {
            return &values_[i * dims_[1]];
        }
        [[nodiscard]] const T* data_ptr_1d(int i) const requires (Dim == 2) {
            return &values_[i * dims_[1]];
        }

        // pointer to column (fixed z) – non-contiguous, stride = Nz

        std::vector<T*> data_ptr_1d_col(size_t j) {
            std::vector<T*> col_ptr(dims_[0]); // Nr elements
            for (size_t i = 0; i < dims_[0]; ++i) {
                col_ptr[i] = &values_[i*dims_[1] + j];
            }
            return col_ptr; // vector of pointers, each pointer points to one element in column z
        }
        // const type version
        [[nodiscard]] std::vector<const T*> data_ptr_1d_col(size_t j) const {
            std::vector<const T*> col_ptr(dims_[0]);
            for (size_t i = 0; i < dims_[0]; ++i)
                col_ptr[i] = &values_[i*dims_[1] + j];
            return col_ptr;
        }

        // raw data access
        std::vector<T> &data() noexcept { return values_; }
        const std::vector<T> &data() const noexcept { return values_; }

        //
        // fill methods
        //

        // fill with index-based callable 2d version
        template<std::invocable<size_t, size_t> F>
        void fill_indexed(F &&f) requires (Dim == 2) {
#pragma omp parallel for collapse(2) default(none) shared(f)
            for (size_t i=0; i<dims_[0]; ++i)
                for (size_t j=0; j<dims_[1]; ++j)
                    values_[linear_idx({i, j})] = std::invoke(f, i, j);
        }

        // fill with index-based callable 1d version
        template<std::invocable<size_t> F>
        void fill_indexed(F &&f) requires (Dim == 1) {
#pragma omp parallel for default(none) shared(f)
            for (size_t i = 0; i < dims_[0]; ++i)
                values_[linear_idx({i})] = std::invoke(f, i);
        }

        // fill with node-based callable 2d (used for r-z slices)
        template<std::invocable<Real, Real> F>
        void fill_nodes(F &&f, std::span<const Real> x_nodes, std::span<const Real> y_nodes) requires (Dim == 2) {
            if (x_nodes.size() != size_t(dims_[0]) || y_nodes.size() != size_t(dims_[1]))
                throw std::invalid_argument("Node vector size mismatch");
#pragma omp parallel for collapse(2) default(none) shared(f, x_nodes, y_nodes)
            for (size_t i = 0; i < dims_[0]; ++i)
                for (size_t j = 0; j < dims_[1]; ++j)
                    values_[linear_idx({i, j})] = std::invoke(f, x_nodes[i], y_nodes[j]);
        }

        // fill with slice with node-based callable 1d (used for z-slices)
        template<std::invocable<Real> F>
        void fill_nodes(F &&f, std::span<const Real> nodes) requires (Dim == 1) {
            if (nodes.size() != size_t(dims_[0]))
                throw std::invalid_argument("Node vector size mismatch");
#pragma omp parallel for default(none) shared(f, nodes)
            for (size_t i = 0; i < dims_[0]; ++i)
                values_[linear_idx({i})] = std::invoke(f, nodes[i]);
        }
    }; // FieldVectorized

/**
 * @brief 1D vectorized GHP field along the "z" direction.
 *
 * @class HeldVectorized
 * HeldFieldVectorized stores an array of Held-scalars with spin weights (p,q) along\n
 * a single row. Supports element-wise arithmetic, conjugation, and
 * node-based filling. Flat contiguous storage enables cache-friendly access.
 * @params
 * Nz   - number of z-points\n
 * init - initial value to fill the field with\n
 * p,q - GHP weights\n
 * @note OpenMP is used to parallelize operations over all elements:\n
 *      - fill_indexed()\n
 *      - fill_nodes()\n
 *      - operator+, operator*\n
 *      - conj()\n
 *
 */
    class HeldVectorized : public FieldVectorized<GHPScalar<Complex>,1> {
    protected:
        int p_{};
        int q_{};
        size_t Nz_{};
    public:
        using Base = FieldVectorized<GHPScalar<Complex>,1>;

        HeldVectorized() = default;
        HeldVectorized(size_t Nz,
                       GHPScalar<Complex> init=GHPScalar<Complex>(teuk::zeroC,0,0),
                       int p=0,int q=0)
                : Base({Nz}, init), p_(p), q_(q), Nz_(Nz){}

        ~HeldVectorized() = default;

        [[nodiscard]] int p() const noexcept { return p_; }
        [[nodiscard]] int q() const noexcept { return q_; }
        [[nodiscard]] size_t Nz() const noexcept { return Nz_; }

        // custom 1d access operator
        GHPScalar<Complex>& operator()(size_t j) { return Base::operator()({j}); }
        const GHPScalar<Complex>& operator()(size_t j) const { return Base::operator()({j}); }

        //
        // GHP preserving arithmetic operators
        //
        HeldVectorized operator+(const HeldVectorized& other) const {
            if (p_ != other.p_ || q_ != other.q_ || dims_[0] != other.dims_[0])
                throw std::runtime_error("Dimension or (p,q) weight mismatch. Cannot add.");
            HeldVectorized out(dims_[0], GHPScalar<Complex>(teuk::zeroC,p_,q_), p_,q_);
#pragma omp parallel for default(none) shared(other, out)
            for (size_t i=0;i<dims_[0];++i) {out(i) = (*this)(i) + other(i);}
            return out;
        }

        HeldVectorized operator*(const HeldVectorized& other) const {
            if (dims_[0] != other.dims_[0]) throw std::runtime_error("Dimension mismatch, cannot multiply.");
            // add GHP weights
            HeldVectorized out(dims_[0], GHPScalar<Complex>(teuk::zeroC,0,0),
                               p_+other.p_, q_+other.q_);
#pragma omp parallel for default(none) shared(other, out)
            for (size_t i=0;i<dims_[0];++i)  {out(i) = (*this)(i) * other(i);}
            return out;
        }

        void conj_inplace() {
            int p0 = p();
            int q0 = q();
            // swap GHP weights
#pragma omp parallel for default(none)
            for (size_t i = 0; i < dims_[0]; ++i)
                (*this)(i).conj();
            assert(p0 == this->q_ && q0 == this->p_ && "weights not swapped correctly");
        }
    }; // HeldFieldVectorized

/**
 * @brief 2D vectorized GHP field container. \n
 * @class GHPVectorized \n
 *
 * GHPFieldVectorized stores a 2D array of GHP-scalars with spin weights
 * (p,q) in a contiguous vector. \n
 * Supports: \n
 *   - row-wise slicing (RSlice) via std::span \n
 *   - column-wise slicing (ZSlice) via pointer arrays \n
 *   - element-wise arithmetic and conjugation \n
 *   - spin-boost transformations \n
 *   - OpenMP-parallelized fill operations \n
 *
 *   @note OpenMP is used to parallelize loops over the 2D field for:  \n
 *    - element-wise operations (operator+, operator*) \n
 *    - filling the field from indexed or nodal functions (fill_indexed, fill_nodes)\n
 *    - conjugation\n
 *    @params
 *    Nr,Nz - number of radial and angular points\n
 *    init   - initial value to fill the field with
 *    p,q   - GHP weights
 *    @details
 *     - ZSlice provides non-contiguous views along fixed z (columns)
 *     - RSlice provides contiguous views along fixed r (rows)
 *
 */
    class GHPVectorized : public FieldVectorized<GHPScalar<Complex>,2> {
    protected:
        int p_ = 0;
        int q_ = 0;
    public:
        using Base = FieldVectorized<GHPScalar<Complex>,2>;
        using Scalar = GHPScalar<Complex>;

        GHPVectorized() = default;
        GHPVectorized(size_t Nr, size_t Nz,
                      GHPScalar<Complex> init=GHPScalar<Complex>(teuk::zeroC, 0, 0),
                      int p=0, int q=0)
                : Base({Nr,Nz}, init), p_(p), q_(q) {}

        [[nodiscard]] int p() const noexcept { return p_; }
        [[nodiscard]] int q() const noexcept { return q_; }


        GHPScalar<Complex>& operator()(size_t r,size_t z) { return Base::operator()({r,z}); }  // operator()(array<int,2>)
        const GHPScalar<Complex>& operator()(size_t r, size_t z) const { return Base::operator()({r,z}); }

        // Apply a unary functor in-place on a span row-slice
        //  span_slice_r(r) gives std::span<GHPScalar<Complex>>
        template<std::invocable<GHPScalar<Complex>&> F>
        void apply_in_place(std::span<GHPScalar<Complex>> sl, F&& f)
        {
#pragma omp parallel for default(none) shared(sl, f)
            for (size_t i = 0; i < sl.size(); ++i)
                std::invoke(f, sl[i]);   // modifies in place
        }

        // ------------------------------------------------------------
        // RSlice structure : views along fixed r
        // ------------------------------------------------------------
        struct RSlice {
            Scalar* data_ptr;   // pointer to start of z-array for fixed r
            size_t Nz_;
            int p_, q_;
            size_t r_;

            RSlice(Scalar* ptr = nullptr, size_t Nz = 0, int p = 0, int q = 0, int r = 0)
                    : data_ptr(ptr), Nz_(Nz), p_(p), q_(q), r_(r) {}

            [[nodiscard]] size_t size() const { return Nz_; }
            [[nodiscard]] int p() const { return p_; }
            [[nodiscard]] int q() const { return q_; }
            [[nodiscard]] size_t r_index() const { return r_; }

            Scalar& operator[](size_t i) {
                assert(data_ptr && i < Nz_);
                return data_ptr[i];
            }
            const Scalar& operator[](size_t i) const {
                assert(data_ptr && i < Nz_);
                return data_ptr[i];
            }
            // conjugate the entire slice (also works with GHP weights by inheriting from GHPScalar)
            void conj_inplace() const {
                for (size_t i = 0; i < Nz_; ++i)
                    data_ptr[i].conj();
            }

        }; // RSlice
// ------------------------------------------------------------
// ZSlice structure : views along fixed z
// ------------------------------------------------------------
        struct ZSlice {
            Scalar* data_ptr; // points to first element at r=0
            size_t Nr_;
            size_t stride_;               // how many elements to jump in memory to next r

            ZSlice(GHPScalar<Complex>* ptr = nullptr, size_t Nr = 0, size_t stride = 0)
                    : data_ptr(ptr), Nr_(Nr), stride_(stride) {}

            [[nodiscard]] size_t size() const { return Nr_; }

            GHPScalar<Complex>& operator[](size_t i) {
                assert(i < Nr_);
                return data_ptr[i*stride_];
            }
            const GHPScalar<Complex>& operator[](size_t i) const {
                assert(i < Nr_);
                return data_ptr[i*stride_];
            }
            void conj_inplace() {
                int p0 = this->data_ptr->p();
                int q0 = this->data_ptr->q();

                for (size_t i=0; i < Nr_; ++i)
                    this->operator[](i).conj();
                // check we changed the weights
                assert(p0==this->data_ptr->q() && q0==this->data_ptr->p());
            }
        };

        // ------------------------------------------------------------
        // Slice accessors
        // ------------------------------------------------------------
        RSlice slice_r(size_t r) {
            Scalar* ptr = Base::data_ptr_1d(r);     // pointer to first element in this row
            return RSlice(ptr, Base::dims_[1], p_, q_, r);
        }
        auto span_slice_r(int r) {
            return std::span<GHPScalar<Complex>>(Base::data_ptr_1d(r), Base::dims_[1]);
        }

        // arithmetic operators
        GHPVectorized operator+(const GHPVectorized& other) const {
            if (p_ != other.p_ || q_ != other.q_ || dims_ != other.dims_)
                throw std::runtime_error("Dimension or (p,q) weight mismatch, cannot add.");
            GHPVectorized out(dims_[0], dims_[1],
                              GHPScalar<Complex>(teuk::zeroC,p_,q_), p_,q_);
#pragma omp parallel for collapse(2) default(none) shared(other, out)
            for (int i=0;i<dims_[0];++i)
                for (int j=0;j<dims_[1];++j)
                    out(i,j) = (*this)(i,j) + other(i,j);
            return out;
        }

        GHPVectorized operator*(const GHPVectorized& other) const {
            if (dims_ != other.dims_) throw std::runtime_error("Dimension mismatch, cannot multiply.");
            GHPVectorized out(dims_[0], dims_[1],
                              GHPScalar<Complex>(teuk::zeroC,0,0),
                              p_+other.p_, q_+other.q_);
#pragma omp parallel for collapse(2) default(none) shared(other, out)
            for (size_t i=0;i<dims_[0];++i)
                for (size_t j=0;j<dims_[1];++j)
                    out(i,j) = (*this)(i,j) * other(i,j);
            return out;
        }

        void conj_inplace() {
#pragma omp parallel for collapse(2) default(none)
            for (size_t i=0;i<dims_[0];++i)
                for (size_t j=0;j<dims_[1];++j)
                    (*this)(i,j) = (*this)(i,j).conj();
        }
    };  // GHPFieldVectorized

/**
 * @name SpectralFieldVectorized
 * @param Nr - number of radial points
 * @param Nz - number of polar points
 * @param modes - struct containing mode indices (m,kr,kz)
 * @param init - initial value to fill the field with
 * @brief 2D spectral field with mode metadata (m,kr,kz).
 *
 * SpectralFieldVectorized is a generic 2D field container for type T \n
 * (e.g., double, complex, GHPScalar ...), augmented with mode metadata: \n
 *   - n: radial mode index\n
 *   - m: azimuthal mode index\n
 *   - k: polar mode index\n
 *
 * Supports flat storage, slicing, resizing, and parallelized fill\n
 * operations.
 */
    template<typename T>
    class SpectralFieldVectorized : public FieldVectorized<T, 2> {

    public:
        struct Modes { int m{}; int kr{}; int kz{}; };

    private:
        bool is_omega_set_ = false;
        bool is_src_set_   = false;
        size_t Nr_{}, Nz_{};
        orbit::SourceFrequencies src_{};
    protected:
        Real omega_mk_{};   // omega_{mk} = m \Omega_\phi  + kr\Omega_r + kz \Omega_z
        Modes modes_;
        // If a derived class *must* mutate modes_ directly, call this after
        void recompute_omega_from_modes() noexcept { update_omega_(); }
    private:
        void update_omega_() noexcept {
            if (!is_src_set_) { is_omega_set_ = false; return; }
            omega_mk_ =
                    Real(modes_.m) *src_.Omega_phi +
                    Real(modes_.kr) * src_.Omega_r   +
                    Real(modes_.kz) * src_.Omega_z;
            is_omega_set_ = true;
        }

    public:
        using Base = FieldVectorized<T, 2>; // build in row-major vectorization functionality
        SpectralFieldVectorized() = default;

        SpectralFieldVectorized(size_t Nr = 0, size_t Nz = 0,
                                Modes modes={0,0,0},
                                const T &init = T())
                : Base({Nr, Nz}, init), modes_(modes), Nr_(Nr), Nz_(Nz) {}
        ~SpectralFieldVectorized() = default;

        // ---- getters ----
        [[nodiscard]] virtual int kr() const noexcept { return modes_.kr; }
        [[nodiscard]] virtual int kz() const noexcept { return modes_.kz; }
        [[nodiscard]] virtual int m() const noexcept { return modes_.m; }
        [[nodiscard]] virtual size_t Nr() const noexcept { return Nr_; }
        [[nodiscard]] virtual size_t Nz() const noexcept { return Nz_; }

        // Prefer using this instead of writing modes_ directly in derived classes.
        void set_modes(Modes md) noexcept {
            modes_ = md;
            update_omega_();
        }
        void set_mode_indices(int m, int kr, int kz) noexcept {
            modes_.m  = m;
            modes_.kr = kr;
            modes_.kz = kz;
            update_omega_();
        }

        [[nodiscard]] const Modes& modes() const noexcept { return modes_; }
        void set_source_frequencies(const orbit::SourceFrequencies& f) noexcept {
            src_ = f;
            is_src_set_ = true;
            update_omega_();
        }

        [[nodiscard]] bool has_source_frequencies() const noexcept { return is_src_set_; }
        [[nodiscard]] const orbit::SourceFrequencies& source_frequencies() const noexcept {
            assert(is_src_set_);
            return src_;
        }
        [[nodiscard]] bool has_omega_mk() const noexcept { return is_omega_set_; }
        [[nodiscard]] Real omega_mk() const noexcept {
            assert(is_omega_set_);
            return omega_mk_;
        }

        [[nodiscard]] Complex iomega_mk() const noexcept {
            assert(is_omega_set_);
            return Complex(0.0, omega_mk_);
        }
        void set_omega_mk(Real omega) noexcept {
            omega_mk_ = omega;
            is_omega_set_ = true;
        }
        void clear_omega_mk() noexcept {
            omega_mk_ = Real(0);
            is_omega_set_ = false;
        }

        virtual T &operator()(size_t ir, size_t iz) { return Base::operator()({ir, iz}); }
        virtual const T &operator()(size_t ir, size_t iz) const { return Base::operator()({ir, iz}); }


        /**
         * @name RSlice
         * @brief non-memory owning view structures for slices\n
         * simply points into storage owned by the parent SpectralFieldVectorized\n
         * mutation operations modify the parent data but no memory is allocated here\n
         * @param data_ptr - pointer to start of slice data (e.g., row for RSlice, column for ZSlice)\n
         * @param Nz_ - number of elements in the slice (e.g., number of z points for RSlice)\n
         * @param r_ - radial index for RSlice (or z index for ZSlice)\n
         * @param modes_ - mode metadata (m, kr, kz) for the slice\n
         * @param omega_mk_ - precomputed frequency combo for the slice (if available)
         *
         */


        struct RSlice {
            T* data_ptr;   // pointer to start of z-array for fixed r
            size_t Nz_;
            size_t r_;
            Modes modes_;

            // per-mode frequency combo omega_mk = m*Omega_phi + kr*Omega_r + kz*Omega_z
            Real omega_mk_{0};
            bool has_omega_{false};

            explicit RSlice(T* ptr=nullptr,
                   size_t Nz=0,
                   Modes modes={0,0,0},
                   size_t r=0, // radial index ***not radial coordinate value ***
                   Real omega=0,
                   bool has_om=false)
                    : data_ptr(ptr), Nz_(Nz), r_(r), modes_(modes),
                      omega_mk_(omega), has_omega_(has_om) {}

            size_t size() const noexcept { return Nz_; }
            size_t r_index()    const noexcept { return r_; }

            [[nodiscard]] int m()  const noexcept { return modes_.m; }
            [[nodiscard]] int kr() const noexcept { return modes_.kr; }
            [[nodiscard]] int kz() const noexcept { return modes_.kz; }

            [[nodiscard]] bool has_omega_mk() const noexcept { return has_omega_; }
            [[nodiscard]] Real omega_mk() const noexcept {
                assert(has_omega_);
                return omega_mk_;
            }
            [[nodiscard]] Complex iomega_mk() const noexcept {
                assert(has_omega_);
                return Complex(0.0, omega_mk_);
            }

            T& operator[](size_t i) {
                assert(data_ptr && i < Nz_);
                return data_ptr[i];
            }
            const T& operator[](size_t i) const {
                assert(data_ptr && i < Nz_);
                return data_ptr[i];
            }

            void conj_inplace() {
                for (size_t i = 0; i < Nz_; ++i) {
                    data_ptr[i] = data_ptr[i].conj();
                }
            }
        }; // RSlice

        // const (view)
        struct ConstRSlice {
            const T* data_ptr;
            size_t Nz_{};
            size_t r_{};
            Modes modes_;

            Real omega_mk_{0};
            bool has_omega_{false};

            const T& operator[](size_t i) const {
                assert(data_ptr && i < Nz_);
                return data_ptr[i];
            }
            size_t size() const noexcept { return Nz_; }
            size_t r_index()    const noexcept { return r_; }

            [[nodiscard]] int m()  const noexcept { return modes_.m; }
            [[nodiscard]] int kr() const noexcept { return modes_.kr; }
            [[nodiscard]] int kz() const noexcept { return modes_.kz; }

            [[nodiscard]] bool has_omega_mk() const noexcept { return has_omega_; }
            [[nodiscard]] Real omega_mk() const noexcept {
                assert(has_omega_);
                return omega_mk_;
            }
            [[nodiscard]] Complex iomega_mk() const noexcept {
                assert(has_omega_);
                return Complex(0.0, omega_mk_);
            }
        };

        struct ZSlice {
            T* data_ptr;            // points to first element at r=0 for fixed z
            size_t Nr_;             // number of r elements on the slice z
            size_t stride_;         // jump to next r (row-major: stride = Nz)
            Modes modes_;
            size_t iz_;

            Real omega_mk_{0};
            bool has_omega_{false};

            ZSlice(T* ptr = nullptr,
                   size_t Nr = 0,
                   size_t stride = 0,
                   Modes modes={0,0,0},
                   size_t iz = 0,
                   Real omega=0,
                   bool has_om=false)
                    : data_ptr(ptr), Nr_(Nr), stride_(stride), modes_(modes), iz_(iz),
                      omega_mk_(omega), has_omega_(has_om) {}

            [[nodiscard]] size_t size() const noexcept { return Nr_; }
            [[nodiscard]] size_t z_index()    const noexcept { return iz_; }

            [[nodiscard]] int kr() const noexcept { return modes_.kr; }
            [[nodiscard]] int kz() const noexcept { return modes_.kz; }
            [[nodiscard]] int m()  const noexcept { return modes_.m; }

            [[nodiscard]] bool has_omega_mk() const noexcept { return has_omega_; }
            [[nodiscard]] Real omega_mk() const noexcept {
                assert(has_omega_);
                return omega_mk_;
            }
            [[nodiscard]] Complex iomega_mk() const noexcept {
                assert(has_omega_);
                return Complex(0.0, omega_mk_);
            }

            T& operator[](size_t i) {
                assert(data_ptr && i < Nr_);
                return data_ptr[i*stride_];
            }
            const T& operator[](size_t i) const {
                assert(data_ptr && i < Nr_);
                return data_ptr[i*stride_];
            }
        }; // ZSlice

        // const (view)
        struct ConstZSlice {
            const T* data_ptr;   // points to first element at r=0 for fixed z
            size_t Nr_;
            size_t stride_;
            Modes modes_;
            size_t iz_;

            // NEW
            Real omega_mk_{0};
            bool has_omega_{false};

            [[nodiscard]] size_t size() const noexcept { return Nr_; }
            [[nodiscard]] size_t z_index()    const noexcept { return iz_; }

            [[nodiscard]] int kr() const noexcept { return modes_.kr; }
            [[nodiscard]] int kz() const noexcept { return modes_.kz; }
            [[nodiscard]] int m()  const noexcept { return modes_.m; }

            [[nodiscard]] bool has_omega_mk() const noexcept { return has_omega_; }
            [[nodiscard]] Real omega_mk() const noexcept {
                assert(has_omega_);
                return omega_mk_;
            }
            [[nodiscard]] Complex iomega_mk() const noexcept {
                assert(has_omega_);
                return Complex(0.0, omega_mk_);
            }

            const T& operator[](size_t i) const {
                assert(data_ptr && i < Nr_);
                return data_ptr[i*stride_];
            }
        };


// =======================
// slice accessors (views)
// =======================

        auto slice_R_span(size_t ir) {
            return std::span<T>(Base::data_ptr_1d(ir), Base::dims()[1]);
        }

        RSlice slice_R(size_t ir) {
            return RSlice(Base::data_ptr_1d(ir),
                          Base::dims()[1],
                          modes_,
                          ir,
                          omega_mk_,
                          is_omega_set_);
        }

        [[nodiscard]] ConstRSlice slice_R(size_t ir) const {
            const T* row = Base::data_ptr_1d(ir);
            return ConstRSlice{ row,
                                Base::dims()[1],
                                ir,
                                modes_,
                                omega_mk_,
                                is_omega_set_ };
        }

// constant z slicing accessors
        //auto slice_Z_span(size_t iz) {
        //    return std::span<T>(Base::data_ptr_1d_col(iz), Base::dims()[0]);
       // }

        ZSlice slice_Z(size_t iz) {
            return ZSlice(Base::values_.data() + iz, // old Base::data_ptr_1d_col(iz).data()+iz, // pointer to first element in column iz
                          Base::dims()[0],
                          Base::dims()[1],
                          Base::modes_,
                          iz,
                          Base::has_omega_mk() ? Base::omega_mk() : Real(0),
                          Base::has_omega_mk());
        }

        ConstZSlice slice_Zconst(size_t iz) const {
            return ConstZSlice{Base::values_.data() + iz, //Base::data_ptr_1d_col(iz)+iz,
                               Base::dims()[0],
                               Base::dims()[1],
                               Base::modes_,
                               iz,
                               Base::has_omega_mk() ? Base::omega_mk() : Real(0),
                               Base::has_omega_mk()};
        }


    }; // SpectralFieldVectorized

/**
 * @brief 2D spectral dim=(Nr)x(Nz) GHP field with spin weights and mode metadata.
 * @class SpectralGHPVectorized
 *
 * @details SpectralGHPVectorized combines the features of GHPFieldVectorized
 * and SpectralFieldVectorized: \n
 *   - stores GHP-scalars with spin weights (p,q) \n
 *   - adds spectral mode metadata (w,m) \n
 *   - supports openMP element-wise arithmetic and conjugation \n
 *   - provides row/column slicing via std::span and pointers \n
 *   - fully contiguous memory layout for cache efficiency \n
 */
    class SpectralGHPVectorized : public SpectralFieldVectorized<GHPScalar<Complex>> {
    protected:
        int p_{};
        int q_{};
    public:
        using Base = SpectralFieldVectorized<GHPScalar<Complex>>;

        SpectralGHPVectorized(size_t Nr, size_t Nz,
                              Modes modes,
                              GHPScalar<Complex> init = GHPScalar<Complex>(
                                      teuk::zeroC, 0, 0), int p=0, int q=0)
                : Base(Nr, Nz, modes, init), p_(p), q_(q)  {}

        SpectralGHPVectorized(size_t Nr, size_t Nz,
                              Modes modes,
                              const ghp::GHPFieldVectorized& field)
                : Base(Nr, Nz, modes, GHPScalar<Complex>(teuk::zeroC, field.p(), field.q())),
                  p_(field.p()), q_(field.q())
        {
            if (field.Nr() != Nr || field.Nz() != Nz) {
                throw std::runtime_error("SpectralGHPVectorized: dimension mismatch");
            }

            for (size_t ir = 0; ir < Nr; ++ir)
                for (size_t iz = 0; iz < Nz; ++iz)
                    (*this)(ir, iz) = field(ir, iz);
        }

        [[nodiscard]] int p() const noexcept { return p_; }
        [[nodiscard]] int q() const noexcept { return q_; }
        [[nodiscard]] int m()  const noexcept { return Base::m(); }
        [[nodiscard]] int kz() const noexcept { return Base::kz(); }
        [[nodiscard]] int kr() const noexcept { return Base::kr(); }
        [[nodiscard]] size_t Nr() const noexcept { return Base::Nr(); }
        [[nodiscard]] size_t Nz() const noexcept { return Base::Nz(); }


        GHPScalar<Complex> &operator()(size_t r, size_t z) { return Base::operator()(r, z); }
        const GHPScalar<Complex> &operator()(size_t r, size_t z) const { return Base::operator()(r, z); }

        inline void set_index(const size_t& ir, const size_t& iz, const GHPScalar<Complex> &val) {
            Base::operator()(ir, iz) = val;
        }

        auto span_slice_r(size_t ir) {
            return std::span<GHPScalar<Complex>>(Base::data_ptr_1d(ir), Base::dims()[1]);
        }

        //Base::RSlice slice_R(size_t r) { return Base::RSlice(Base::data_ptr_1d(r), Base::dims()[1],
         //                                                    Base::modes_, r); }

        Base::RSlice slice_R(size_t r) {
            return Base::RSlice(Base::data_ptr_1d(r),
                                Base::dims()[1],
                                Base::modes_,
                                r,
                                Base::has_omega_mk() ? Base::omega_mk() : Real(0),
                                Base::has_omega_mk());
        }
        Base::ConstRSlice slice_R(size_t r) const {
            return Base::ConstRSlice{
                    Base::data_ptr_1d(r),
                    Base::dims()[1],
                    r,
                    Base::modes_,
                    Base::has_omega_mk() ? Base::omega_mk() : Real(0),
                    Base::has_omega_mk()
            };
        }

        std::vector<GHPScalar<Complex>> col_slice_z(size_t iz) const {
            std::vector<GHPScalar<Complex>> col(Base::dims()[0]);
            for (size_t r = 0; r < Base::dims()[0]; ++r) {
                col[r] = (*this)(r, iz);
            }
            return col; // copy of the column
        }

        //Base::ZSlice slice_Z(size_t iz) {
        //    return Base::ZSlice(Base::values_.data() + iz, Base::dims()[0], Base::dims()[1], Base::modes_, iz);
        //}
        Base::ZSlice slice_Z(size_t iz) {
            return Base::ZSlice(Base::values_.data() + iz,
                                Base::dims()[0],
                                Base::dims()[1],
                                Base::modes_,
                                iz,
                                Base::has_omega_mk() ? Base::omega_mk() : Real(0),
                                Base::has_omega_mk());
        }
        Base::ConstZSlice slice_Z(size_t iz) const {
            return Base::ConstZSlice(Base::values_.data() + iz,
                                Base::dims()[0],
                                Base::dims()[1],
                                Base::modes_,
                                iz,
                                Base::has_omega_mk() ? Base::omega_mk() : Real(0),
                                Base::has_omega_mk());
        }
        // arithmetic operators
        SpectralGHPVectorized operator+(const SpectralGHPVectorized &other) const {
            if (p_ != other.p_ || q_ != other.q_ ||
                Base::dims() != other.Base::dims() || Base::m() != other.m() ||
                Base::kr() != other.kr() || Base::kz() != other.kz())
                throw std::runtime_error("Mismatch in weights, dims or modes");
            SpectralGHPVectorized out(Base::dims()[0], Base::dims()[1], Base::modes_,
                                      GHPScalar<Complex>(teuk::zeroC, p_, q_), p_, q_);
#pragma omp parallel for collapse(2) default(none) shared(other, out)
            for (size_t i = 0; i < Base::dims()[0]; ++i)
                for (size_t j = 0; j < Base::dims()[1]; ++j)
                    out(i, j) = (*this)(i, j) + other(i, j);
            return out;
        }

        SpectralGHPVectorized operator*(const SpectralGHPVectorized &other) const {
            if (Base::dims() != other.Base::dims()) throw std::runtime_error("Dimension mismatch");
            SpectralGHPVectorized out(Base::dims()[0], Base::dims()[1],  Base::modes_,
                                      GHPScalar<Complex>(teuk::zeroC, 0, 0), p_ + other.p_, q_ + other.q_);
#pragma omp parallel for collapse(2) default(none) shared(other, out)
            for (size_t i = 0; i < Base::dims()[0]; ++i)
                for (size_t j = 0; j < Base::dims()[1]; ++j)
                    out(i, j) = (*this)(i, j) * other(i, j);
            return out;
        }

        void conj_inplace() {
#pragma omp parallel for collapse(2) default(none)
            for (size_t i = 0; i < Base::dims()[0]; ++i)
                for (size_t j = 0; j < Base::dims()[1]; ++j)
                    (*this)(i, j).conj();
        }

    }; // SpectralGHPVectorized

}; // spectral

#endif // GHZ_NUMERIC_SPECTRALGHPFIELDVECTORIZED_HPP