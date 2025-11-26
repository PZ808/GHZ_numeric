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
// Generic container storing data_[omega][m][r][z] where each scalar is type T
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

        const Scalar &at(int r, int z) const { return values_(r, z); }

        //
        // ZSlice: non-owning view along z for fixed r (for this GHPSpectralField)
        //
        struct ZSlice {
            Scalar *data_ptr; // pointer to first z element at this r
            int Nz_;
            int w_, m_, r_;

            ZSlice(Scalar *ptr = nullptr, int Nz = 0, int w = 0, int m = 0, int r = 0)
                    : data_ptr(ptr), Nz_(Nz), w_(w), m_(m), r_(r) {}

            int size() const { return Nz_; }
            int m() const { return m_; }
            int w() const { return w_; }
            int r() const { return r_; }

            Scalar &operator[](int i) {
                assert(data_ptr && i >= 0 && i < Nz_);
                return data_ptr[i];
            }

            const Scalar &operator[](int i) const {
                assert(data_ptr && i >= 0 && i < Nz_);
                return data_ptr[i];
            }
        }; // ZSlice

        // Return a non-owning ZSlice for fixed r
        ZSlice slice_z(int r) {
            auto &row = values_.row_r(r);
            return ZSlice(row.data(), values_.Nz(), w_, m_, r);
        }

        // const version
        ZSlice slice_z(int r) const {
            auto &row = values_.row_r(r);
            return ZSlice(const_cast<Scalar *>(row.data()), values_.Nz(), w_, m_, r);
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
            GHPSpectralField out(w_, m_, Nr(), Nz(), Scalar(teuk::zeroC, p_, q_), p_, q_);
            for (int r = 0; r < Nr(); ++r)
                for (int z = 0; z < Nz(); ++z)
                    out.at(r, z) = this->at(r, z) + other.at(r, z);
            return out;
        }

        GHPSpectralField transform(const Complex &lambda) const {
            GHPSpectralField out(w_, m_, Nr(), Nz(), Scalar(teuk::zeroC, p_, q_), p_, q_);
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



// -----------------------------------------------------------------------------
// SpectralDiffer: operates on Z-slices for a given block (fixed w,m)
// -----------------------------------------------------------------------------
    class SpectralDiffer {
    private:
        int Nz_;                       // Number of LGL nodes (theta)
        bool use_fast_double_fft_;      // Whether to use double FFT

        std::vector<Real> z_;           // LGL nodes [-1,1]
        std::vector<Real> w_;           // LGL quadrature weights
        matrix<Real> D_;                // Legendre differentiation matrix

    public:
        // Constructor: build nodes, weights, differentiation matrix
        SpectralDiffer(int Nz, bool use_fast_double_fft = true)
                : Nz_(Nz), use_fast_double_fft_(use_fast_double_fft) {
            z_ = build_legendre_gauss_lobatto_nodes(Nz_); // generate nodes z[i]
            w_ = build_barycentric_weights();
            D_ = build_legendre_diff_matrix(z_);
        }

        // Getters
        const std::vector<Real> &nodes() const { return z_; }
        const std::vector<Real> &weights() const { return w_; }
        const matrix<Real> &Dmatrix() const { return D_; }

        // Barycentric utilities (declare; implement elsewhere or below)
        std::vector<Real> build_barycentric_weights() const;
        std::pair<Complex, Complex> barycentric_interp_and_derivative(
                const std::vector<Complex> &f, const std::vector<Real> &w, Real z0) const;

        static std::pair<Real, Real> legendre_P_and_dP(int n, Real x);
        static std::vector<Real> build_legendre_gauss_lobatto_nodes(int N);
        static matrix<Real> build_legendre_diff_matrix(const std::vector<Real> &x);

        // dz on a single ZSlice (in-place into df_dz)
        void dz_Dmatrix(const GHPSpectralField::ZSlice &f, GHPSpectralField::ZSlice &df_dz) const;
        GHPSpectralField::ZSlice dz_barycentric(const GHPSpectralField::ZSlice &f, GHPSpectralField::ZSlice &df_dz,
                                                const std::vector<Real> &w) const;

        // d/dphi for a single slice
        void dphi_fft_single_m(const GHPSpectralField::ZSlice &f_slice,
                               GHPSpectralField::ZSlice &df_dphi_slice) const;

        // edth operator on a single slice -> returns a new ZSlice-owned buffer (or you can provide out)
        GHPSpectralField::ZSlice edth(const GHPSpectralField::ZSlice &f) const;


    };

}
#endif // SPECTRAL_GHPFIELD_HPP