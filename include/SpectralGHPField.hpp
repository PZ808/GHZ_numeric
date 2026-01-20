//
// Created by Peter Zimmerman on 22.11.25.
//

#ifndef GHZ_NUMERIC_SPECTRALGHPFIELD_HPP
#define GHZ_NUMERIC_SPECTRALGHPFIELD_HPP


#include <vector>
#include <complex>
#include <functional>
#include <stdexcept>
#include <sstream>
#include <iomanip>
#include <cassert>
#include <complex>
#include <type_traits>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include "GhzTypes.hpp"
#include "GHPScalars.hpp"
#include "HeldScalars.hpp"
#include "VectorsGHZ.hpp"  // Your fixed-size vector class

namespace spectral {

    using Complex = teuk::Complex;
    using Real = teuk::Real;
    using boost::numeric::ublas::matrix;

// -----------------------------------------------------------------------------
// SpectralField<T>
// Generic container storing data_[r][z]
// -----------------------------------------------------------------------------
    template<typename T>
    class SpectralField {
    public:
        using ValueType = T;
        using ZRow = std::vector<T>;            // z array for a single r
        using RZBlock = std::vector<ZRow>;      // [r][z], where r is a physical grid index and z a spectral index

    private:
        int w_, m_, Nr_, Nz_;
        RZBlock data_;

    public:
        // Constructor
        SpectralField(int Nr = 0, int Nz = 0, int m = 0, int w = 0, const T &init = T())
                : Nr_(Nr), Nz_(Nz), m_(m), w_(w),
                data_(RZBlock(Nr, std::vector<T>(Nz, init))) {}

        // Resize
        void resize(int Nr, int Nz, const T &init = T()) {
            Nr_ = Nr;
            Nz_ = Nz;
            data_.assign(RZBlock(Nr, std::vector<T>(Nz, init)));
        }

        // modes
        int w() const noexcept { return w_; }
        int m() const noexcept { return m_; }
        // Dimensions
        int Nr() const noexcept { return Nr_; }
        int Nz() const noexcept { return Nz_; }

        // access a full r-row [z]
        ZRow &row_r(int r) {
            assert(r >= 0 && r < Nr_);
            return data_[r];
        }

        const ZRow &row_r(int r) const {
            assert(r >= 0 && r < Nr_);
            return data_[r];
        }

        // Element access
        T &operator()(int r, int z) {
            assert(r >= 0 && r < Nr_ && z >= 0 && z < Nz_);
            return data_[r][z];
        }

        const T &operator()(int r, int z) const {
            assert(r >= 0 && r < Nr_ && z >= 0 && z < Nz_);
            return data_[r][z];
        }

        // fill with callable f(r,z) -> T using lambda function
        void fill(std::function<T(int,int)> func) {
            for (int r = 0; r < Nr_; ++r)
                for (int z = 0; z < Nz_; ++z)
                    data_[r][z] = func(r,z);
        }

        // block-wise conjugation for this single block (requires T::conj())
        void conj_block() {
            for (int r = 0; r < Nr_; ++r)
                for (int z = 0; z < Nz_; ++z) {
                    if constexpr (std::is_member_function_pointer<decltype(&T::conj)>::value) {
                        data_[r][z] = data_[r][z].conj();
                    } else {
                        data_[r][z] = std::conj(data_[r][z]);
                    }
                }
        }
    };
/**
* \class GHPSpectralField
* \parameters
* w = omega, m = azimode #,
* Nz = numbers of polar points, Nr = number of radial grid points
* values = array of initial values
*  p = GHP weight, q = GHP weight
*/
    class GHPSpectralField {
    public:
        using Scalar = GHPScalar<Complex>;
        using Field = SpectralField<Scalar>;

    private:
        Field values_; // values_[omega][m][r][z]
        int m_, w_, p_, q_;    // spin weights

    public:
        GHPSpectralField() : values_(), w_(0), m_(0), p_(0), q_(0) {}

        GHPSpectralField(int w, int m, int Nr, int Nz,
                         Scalar init = Scalar(teuk::zeroC, 0, 0),
                         int p = 0, int q = 0)
                : values_(Nr, Nz, m, w, init), p_(p), q_(q) {
            // ensure each scalar carries correct spin weights
            values_.fill([p, q](int r, int z) {
                return Scalar(teuk::zeroC, p, q);
            });
        }

        // metadata & dims
        int w() const noexcept { return w_; }
        int m() const noexcept { return m_; }
        int Nr() const noexcept { return values_.Nr(); }
        int Nz() const noexcept { return values_.Nz(); }
        int p() const noexcept { return p_; }
        int q() const noexcept { return q_; }

        // element access
        Scalar &at(int r, int z) { return values_(r, z); }

        [[nodiscard]] const Scalar &at(int r, int z) const { return values_(r, z); }

        //
        // RSlice: non-owning view along z for fixed r (for this GHPSpectralField)
        //
        //
        struct RSlice {
            Scalar* data_ptr;   // pointer to first z element at this r
            int Nz_;
            int w_, m_, r_;

            RSlice(Scalar* ptr = nullptr,
                   int Nz=0,
                   int w=0,
                   int m=0,
                   int r=0)
                    : data_ptr(ptr), Nz_(Nz), w_(w), m_(m), r_(r)
            {}

            // ---------------------------------------------------------------------
            // Basic info
            // ---------------------------------------------------------------------
            int size() const { return Nz_; }
            int w() const { return w_; }
            int m() const { return m_; }
            int r() const { return r_; }

            // ---------------------------------------------------------------------
            // Element access
            // ---------------------------------------------------------------------
            Scalar& operator[](int i) {
                assert(data_ptr && i >= 0 && i < Nz_);
                return data_ptr[i];
            }

            const Scalar& operator[](int i) const {
                assert(data_ptr && i >= 0 && i < Nz_);
                return data_ptr[i];
            }

            // ---------------------------------------------------------------------
            // Iterator (contiguous memory → true random-access iterator)
            // ---------------------------------------------------------------------
            struct iterator {
                Scalar* ptr;

                using iterator_category = std::random_access_iterator_tag;
                using value_type        = Scalar;
                using reference         = Scalar&;
                using pointer           = Scalar*;
                using difference_type   = std::ptrdiff_t;

                explicit iterator(Scalar* p = nullptr) : ptr(p) {}

                reference operator*() const { return *ptr; }
                pointer operator->() const { return ptr; }

                iterator& operator++() { ptr++; return *this; }
                const iterator operator++(int) { iterator tmp = *this; ptr++; return tmp; }

                iterator& operator--() { ptr--; return *this; }
                const iterator operator--(int) { iterator tmp = *this; ptr--; return tmp; }

                iterator operator+(difference_type n) const { return iterator(ptr + n); }
                iterator operator-(difference_type n) const { return iterator(ptr - n); }
                difference_type operator-(const iterator& other) const { return ptr - other.ptr; }

                bool operator==(const iterator& o) const { return ptr == o.ptr; }
                bool operator!=(const iterator& o) const { return ptr != o.ptr; }
                bool operator<(const iterator& o) const { return ptr < o.ptr; }
                bool operator>(const iterator& o) const { return ptr > o.ptr; }
                bool operator<=(const iterator& o) const { return ptr <= o.ptr; }
                bool operator>=(const iterator& o) const { return ptr >= o.ptr; }
            };

            struct const_iterator {
                const Scalar* ptr;

                using iterator_category = std::random_access_iterator_tag;
                using value_type        = Scalar;
                using reference         = const Scalar&;
                using pointer           = const Scalar*;
                using difference_type   = std::ptrdiff_t;

                explicit const_iterator(const Scalar* p = nullptr) : ptr(p) {}

                reference operator*() const { return *ptr; }
                pointer operator->() const { return ptr; }

                const_iterator& operator++() { ptr++; return *this; }
                const const_iterator operator++(int) { const_iterator tmp = *this; ptr++; return tmp; }

                const_iterator& operator--() { ptr--; return *this; }
                const const_iterator operator--(int) { const_iterator tmp = *this; ptr--; return tmp; }

                const_iterator operator+(difference_type n) const { return const_iterator(ptr + n); }
                const_iterator operator-(difference_type n) const { return const_iterator(ptr - n); }
                difference_type operator-(const const_iterator& other) const { return ptr - other.ptr; }

                bool operator==(const const_iterator& o) const { return ptr == o.ptr; }
                bool operator!=(const const_iterator& o) const { return ptr != o.ptr; }
            };

            // ---------------------------------------------------------------------
            // begin/end
            // ---------------------------------------------------------------------
            iterator begin() { return iterator(data_ptr); }
            iterator end()   { return iterator(data_ptr + Nz_); }

            const_iterator begin() const { return const_iterator(data_ptr); }
            const_iterator end()   const { return const_iterator(data_ptr + Nz_); }

            const_iterator cbegin() const { return const_iterator(data_ptr); }
            const_iterator cend()   const { return const_iterator(data_ptr + Nz_); }

            // ---------------------------------------------------------------------
            // std::span support (contiguous = allowed!)
            // ---------------------------------------------------------------------
            std::span<Scalar> span() {
                return std::span<Scalar>(data_ptr, Nz_);
            }

            std::span<const Scalar> span() const {
                return std::span<const Scalar>(data_ptr, Nz_);
            }

            // ---------------------------------------------------------------------
            // Conjugated copy
            // ---------------------------------------------------------------------
            RSlice conj() const {
                RSlice out(data_ptr, Nz_, w_, m_, r_);
                for (int i = 0; i < Nz_; ++i)
                    out[i] = this->operator[](i).conj();
                return out;
            }
        }; // RSlice

        struct ConstRSlice {
            const Scalar* data_ptr;
            int Nz_, w_, m_, r_;

            const Scalar& operator[](int i) const { return data_ptr[i]; }
            int size() const { return Nz_; }
        };

        ConstRSlice slice_r(int r) const {
            auto &row = values_.row_r(r);
            return ConstRSlice{row.data(), values_.Nz(), w_, m_, r};
        }
        // Return a non-owning RSlice for fixed r
        RSlice slice_r(int r) {
            auto &row = values_.row_r(r);
            return RSlice(row.data(), values_.Nz(), w_, m_, r);
        }

        // const version
        //RSlice slice_r(int r) const {
        //    auto &row = values_.row_r(r);
        //    return RSlice(const_cast<Scalar *>(row.data()), values_.Nz(), w_, m_, r);
        //}

        struct ZSlice {
            Scalar* data_ptr;    // pointer to element at (0, z)
            int Nr_;             // number of radial points
            int Nz_;             // stride in memory (full z dimension)
            int w_, m_, p_, q_;  // GHP / mode metadata
            int z_;              // fixed column index

            ZSlice(Scalar* ptr = nullptr,
                   int Nr = 0, int Nz = 0,
                   int w = 0, int m = 0, int p = 0, int q = 0,
                   int z = 0)
                    : data_ptr(ptr), Nr_(Nr), Nz_(Nz),
                      w_(w), m_(m), p_(p), q_(q), z_(z)
            {}

            // --- size & metadata ------------------------------------------------------
            int size() const { return Nr_; }
            int m() const { return m_; }
            int w() const { return w_; }
            int p() const { return p_; }
            int q() const { return q_; }
            int z() const { return z_; }

            // --- element access -------------------------------------------------------
            Scalar& operator[](int r) {
                assert(data_ptr && r >= 0 && r < Nr_);
                return *(data_ptr + r * Nz_);
            }

            const Scalar& operator[](int r) const {
                assert(data_ptr && r >= 0 && r < Nr_);
                return *(data_ptr + r * Nz_);
            }

            // --- iterator (forward, uses stride) -------------------------------------
            struct iterator {
                Scalar* ptr;
                int stride;

                using iterator_category = std::forward_iterator_tag;
                using value_type        = Scalar;
                using difference_type   = std::ptrdiff_t;
                using pointer           = Scalar*;
                using reference         = Scalar&;

                iterator(Scalar* p, int s) : ptr(p), stride(s) {}

                reference operator*() const { return *ptr; }
                pointer operator->() const { return ptr; }

                iterator& operator++() { ptr += stride; return *this; }
                iterator operator++(int) { iterator tmp = *this; ptr += stride; return tmp; }

                bool operator==(const iterator& other) const { return ptr == other.ptr; }
                bool operator!=(const iterator& other) const { return ptr != other.ptr; }
            };

            struct const_iterator {
                const Scalar* ptr;
                int stride;

                using iterator_category = std::forward_iterator_tag;
                using value_type        = Scalar;
                using difference_type   = std::ptrdiff_t;
                using pointer           = const Scalar*;
                using reference         = const Scalar&;

                const_iterator(const Scalar* p, int s) : ptr(p), stride(s) {}

                reference operator*() const { return *ptr; }
                pointer operator->() const { return ptr; }

                const_iterator& operator++() { ptr += stride; return *this; }
                const_iterator operator++(int) { const_iterator tmp = *this; ptr += stride; return tmp; }

                bool operator==(const const_iterator& other) const { return ptr == other.ptr; }
                bool operator!=(const const_iterator& other) const { return ptr != other.ptr; }
            };

            iterator begin() { return iterator(data_ptr, Nz_); }
            iterator end()   { return iterator(data_ptr + Nr_ * Nz_, Nz_); }

            const_iterator begin() const { return const_iterator(data_ptr, Nz_); }
            const_iterator end()   const { return const_iterator(data_ptr + Nr_ * Nz_, Nz_); }

            const_iterator cbegin() const { return begin(); }
            const_iterator cend()   const { return end(); }

            // --- conj copy ------------------------------------------------------------
            ZSlice conj() const {
                ZSlice out(data_ptr, Nr_, Nz_, w_, m_, p_, q_, z_);
                for (int r = 0; r < Nr_; ++r)
                    out[r] = this->operator[](r).conj();
                return out;
            }

            // === Strided span proxy ==================================================
            //
            // True std::span cannot represent strided memory. This helper gives you the
            // same interface but internally uses the iterator with stride.
            //
            struct Span {
                ZSlice& parent;

                Span(ZSlice& zs) : parent(zs) {}

                int size() const { return parent.size(); }
                Scalar& operator[](int i) { return parent[i]; }
                const Scalar& operator[](int i) const { return parent[i]; }

                iterator begin() { return parent.begin(); }
                iterator end()   { return parent.end(); }
            };

            struct ConstSpan {
                const ZSlice& parent;

                ConstSpan(const ZSlice& zs) : parent(zs) {}

                int size() const { return parent.size(); }
                const Scalar& operator[](int i) const { return parent[i]; }

                const_iterator begin() const { return parent.begin(); }
                const_iterator end()   const { return parent.end(); }
            };

            Span span() { return Span(*this); }
            ConstSpan span() const { return ConstSpan(*this); }

        }; //ZSlice



        ZSlice slice_z(int z) {
            Scalar* ptr = &values_(0, z); // pointer to (0,z)
            return ZSlice(ptr, Nr(), Nz(), w_, m_, p_, q_, z);
        }

        // Element-wise operations (mode-to-mode: require same w,m,Nr,Nz)
        GHPSpectralField conj() const {
            GHPSpectralField out(w_, m_, Nr(), Nz(), Scalar(teuk::zeroC, 0, 0), q_, p_);
            for (int r = 0; r < Nr(); ++r)
                for (int z = 0; z < Nz(); ++z)
                    out.at(r, z) = this->at(r, z).conj();
            return out;
        }

        GHPSpectralField operator*(const GHPSpectralField &other) const {
            if (w_ != other.w_ || m_ != other.m_ || Nr() != other.Nr() || Nz() != other.Nz())
                throw std::runtime_error("GHPSpectralField mismatch in multiplication");
            GHPSpectralField out(w_, m_, Nr(), Nz(), Scalar(teuk::zeroC, 0, 0), p_ + other.p_, q_ + other.q_);
            for (int r = 0; r < Nr(); ++r)
                for (int z = 0; z < Nz(); ++z)
                    out.at(r, z) = this->at(r, z) * other.at(r, z);
            return out;
        }

        // multiply each element by a GHPScalar constant (e.g. i*m as a scalar)
        GHPSpectralField operator*(const GHPScalar<Complex> &scalar) const {
            GHPSpectralField out(w_, m_, Nr(), Nz(), Scalar(teuk::zeroC, 0, 0),
                                 p_ + scalar.p(), q_ + scalar.q());
            for (int r = 0; r < Nr(); ++r)
                for (int z = 0; z < Nz(); ++z)
                    out.at(r, z) = this->at(r, z) * scalar;
            return out;
        }

        GHPSpectralField operator+(const GHPSpectralField &other) const {
            if (w_ != other.w_ || m_ != other.m_ || Nr() != other.Nr() || Nz() != other.Nz())
                throw std::runtime_error("GHPSpectralField mismatch in addition");
            if (p_ != other.p_ || q_ != other.q_)
                throw std::runtime_error("GHPSpectralField spin-weights mismatch in addition");
            GHPSpectralField out(w_, m_,
                                 Nr(), Nz(),
                                 Scalar(teuk::zeroC, p_, q_),
                                 p_, q_);
            for (int r = 0; r < Nr(); ++r)
                for (int z = 0; z < Nz(); ++z)
                    out.at(r, z) = this->at(r, z) + other.at(r, z);
            return out;
        }

        GHPSpectralField transform(const Complex &lambda) const {
            GHPSpectralField out(w_, m_,
                                 Nr(), Nz(),
                                 Scalar(teuk::zeroC, p_, q_),
                                 p_, q_);
            for (int r = 0; r < Nr(); ++r)
                for (int z = 0; z < Nz(); ++z)
                    out.at(r, z) = this->at(r, z).transform(lambda);
            return out;
        }

        // fill from callable f(r,z) -> Complex (wrap to Scalar with p,q)
        void fill(std::function<Complex(int, int)> func) {
            for (int r = 0; r < Nr(); ++r)
                for (int z = 0; z < Nz(); ++z)
                    values_(r, z) = Scalar(func(r, z), p_, q_);
        }

        // debug string
        std::string str(int r = -1, int z = -1) const {
            std::ostringstream oss;
            oss << "GHPSpectralField (w,m)=(" << w_ << "," << m_ << "), (p,q)=(" << p_ << "," << q_
                << "), size=(Nr,Nz)=(" << Nr() << "," << Nz() << ")\n";
            if (r >= 0 && z >= 0) {
                const auto &s = at(r, z);
                oss << "value = " << std::setprecision(8) << s.value()
                    << "  (p,q)=(" << s.p() << "," << s.q() << ")\n";
            }
            return oss.str();
        }
    };

}
#endif // GHZ_NUMERIC_SPECTRALGHPFIELD_HPP