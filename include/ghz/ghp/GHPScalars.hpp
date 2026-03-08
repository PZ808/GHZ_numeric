//
// Created by Peter Zimmerman on 25.10.25.
//

#ifndef GHZ_NUMERIC_GHPSCALARS_HPP
#define GHZ_NUMERIC_GHPSCALARS_HPP

#include "SpinCoeffsNP.hpp"
#include "ghz/core/MathMacros.hpp"
#include "SpinCoeffsNP.hpp"

#include <complex>
#include <iomanip>

#ifdef _OPENMP
#include <omp.h>
#endif
using namespace math;

/**
 * @class GHPScalar
 * @brief Represents a Geroch-Held-Penrose (GHP) scalar with spin-boost weights (p,q)
 *
 * @details This class wraps a complex number (or any ComplexT type) as a scalar in the GHP formalism. \n
 * Each scalar carries spin-boost weights (p,q) which define its transformation properties under \n
 * null tetrad rotations and boosts. \n
 *
 * Features: \n
 * - Access and modify the underlying complex value. \n
 * - Retrieve spin-boost weights and (p,q). \n
 * - GHP conjugation: flips weights (p,q) -> (q,p) and complex conjugates the value. \n
 * - Element-wise arithmetic operations that properly combine weights: \n
 *     - Multiplication: weights add \n
 *     - Division: weights subtract \n
 *     - Addition/Subtraction: only allowed for matching weights \n
 * - Spin-boost transformation: multiplies the value by λ^p * λ̄^q \n
 *
 * @tparam ComplexT Type of the underlying complex number (default: teuk::Complex) used for boost multi
 */
template <typename ComplexT = teuk::Complex>
class GHPScalar {
public:
    using Complex = ComplexT;
    using Real = teuk::Real;

protected:
    int p_;
    int q_;
private:
    Complex value_;
    //int p_; // weight
    //int q_; // weight
    int spin_;
    int boost_;

public:
    explicit GHPScalar(Complex val = teuk::zeroC, int p = 0, int q = 0)
            : value_(val), p_(p), q_(q)
    {
        boost_ = (p_+q_)/2; spin_ = (p_-q_)/2;
        if ( ((p+q) % 2 != 0) || ((p-q) % 2 != 0) ) {
            throw std::runtime_error("Invalid GHP weights: spin/boost would not be integers");
        }
    }

    [[nodiscard]] Complex value() const { return value_; }
    [[nodiscard]] Complex& value() { return value_; }
    [[nodiscard]] int p() const { return p_; }
    [[nodiscard]] int q() const { return q_; }
    [[nodiscard]] int s() const { return spin_; }
    [[nodiscard]] int b() const { return boost_; }

    void set_pq(int p, int q) { p_ = p; q_ = q; boost_ = (int)(p_+q_)/2; spin_ = (int)(p_-q_)/2;}
    void setValue(Complex val) { value_ = val; }

    // conjugation flips weights (p,q) -> (q,p)
    GHPScalar conj() const { return GHPScalar(std::conj(value_), q_, p_); }
    // multiplication combines weights
     GHPScalar operator*(const GHPScalar& other) const {
        return GHPScalar(value_ * other.value_, p_ + other.p_, q_ + other.q_);
    }
    GHPScalar operator/(const GHPScalar& other) const {
        return GHPScalar(value_ / other.value_, p_ - other.p_, q_ - other.q_);
    }
    GHPScalar operator+(const GHPScalar& other) const {
        // addition only makes sense if weights match
        if (p_ != other.p_ || q_ != other.q_) {
            throw std::runtime_error("Incompatible GHP weights in addition");
        }
        return GHPScalar(value_ + other.value_, p_, q_);
    }
    GHPScalar operator-(const GHPScalar& other) const {
        // addition only makes sense if weights match
        if (p_ != other.p_ || q_ != other.q_) {
            throw std::runtime_error("Incompatible GHP weights in addition");
        }
        return GHPScalar(value_ - other.value_, p_, q_);
    }
    // unary minus
    GHPScalar operator-() const {
        return GHPScalar(-value_, p_, q_);
    }
    // in-place operators
    GHPScalar& operator+=(const GHPScalar& other) {
        if (p_ != other.p_ || q_ != other.q_) {
            throw std::runtime_error("Incompatible GHP weights in addition");
        }
        value_ += other.value_;
        return *this;
    }
    GHPScalar& operator-=(const GHPScalar& other) {
        if (p_ != other.p_ || q_ != other.q_) {
            throw std::runtime_error("Incompatible GHP weights in subtraction");
        }
        value_ -= other.value_;
        return *this;
    }
    template<typename U>
    GHPScalar& operator*=(U a) {
        value_ *= a;
        return *this;
    }
    template<typename U>
    GHPScalar& operator/=(U a) {
        value_ /= a;
        return *this;
    }
    bool operator==(const GHPScalar& other) const {
        return p_ == other.p_ && q_ == other.q_ && value_ == other.value_;
    }

    bool operator!=(const GHPScalar& other) const {
        return !(*this == other);
    }

    // multiplying ordinary numbers
    template<typename U>
    requires std::is_arithmetic_v<U> || teuk::is_complex_v<U> || teuk::is_teuk_scalar_v<U>
    auto operator*(U a) const {
        using R = decltype(value_ * a);
        return GHPScalar<R>(value_ * a, p_, q_);
    }
    template<typename U>
    requires std::is_arithmetic_v<U> || teuk::is_teuk_scalar_v<U>
    friend auto operator*(U a, const GHPScalar& x) { return x * a; }

    // dividing by ordinary numbers
    template<typename U>
    requires std::is_arithmetic_v<U> || teuk::is_teuk_scalar_v<U>
    auto operator/(U a) const {
        using R = decltype(value_ / a);
        return GHPScalar<R>(value_ / a, p_, q_);
    }
    // dividing ordinary number by GHP scalar
    template<typename U>
    requires std::is_arithmetic_v<U> || teuk::is_teuk_scalar_v<U>
    friend auto operator/(U a, const GHPScalar& x) {
        using R = decltype(a / x.value_);
        return GHPScalar<R>(a / x.value_, -x.p_, -x.q_);
    }

    // spin-boost transformation
    [[nodiscard]] GHPScalar transform(const Complex& lambda) const {
        Complex newVal = math::PowInt(lambda, p_) * math::PowInt(std::conj(lambda), q_) * value_;
    return GHPScalar(newVal, p_, q_);
    }
    // debug string
    std::string str() const {
        return "val=" + std::to_string(value_.real()) +
               (value_.imag() >= 0 ? "+" : "") + std::to_string(value_.imag()) + "i"
               + " (p,q)=(" + std::to_string(p_) + "," + std::to_string(q_) + ")";
    }
}; // class GHPScalar

/**
 * @class GHPField
 * @brief 2D spectral/storage container for a Geroch–Held–Penrose (GHP) scalar field.
 *
 * This class stores a rectangular grid of GHPScalar values representing a field
 * Φ(r,z) with fixed GHP weights (p,q). It provides:
 *
 *  - element access via operator()(r,z)
 *  - element-wise addition/multiplication of fields
 *  - conjugation (GHP prime operation)
 *  - spin/boost transformation
 *  - initialization and callable-based filling
 *
 * The underlying data layout is:
 *      values_[r][z]   where:
 *          r = 0 ... Nr-1
 *          z = 0 ... Nz-1
 *
 * Weights (p,q) track the GHP type of the field and propagate automatically
 * through algebraic operations.
 *
 * Typical usage:
 *
 *      GHPField phi(Nr, Nz, Scalar(0, p, q), p, q);
 *      phi.fill([](int r,int z){ return complex<double>(r+z,0); });
 *      auto psi = phi.conj();   // GHP conjugation swaps (p,q)
 *
 * This class does not own any coordinate information — it only stores values.
 * Coordinate grids (r,z nodes) are passed to the fill(...) overload if needed.
 */
class GHPField {
public:
    using Complex = teuk::Complex;
    using Scalar = GHPScalar<Complex>;
    using RZBlockScalar = std::vector<std::vector<Scalar>>;

private:
    RZBlockScalar values_; //  Nr x Nz array
    int p_;
    int q_;
    size_t Nz_;
    size_t Nr_;
    Scalar zero_pq_;

public:
    /**
    * @brief Construct an Nr × Nz GHP field with fixed GHP weights (p,q).
    *
    * @param Nr   Number of radial points
    * @param Nz   Number of z (angular/collocation) points
    * @param init Initial scalar used to fill the field (value only; p,q overwritten)
    * @param p    GHP weight p
    * @param q    GHP weight q
    *
    * All entries are initialized as GHPScalar(init.value(), p, q).
    */
    GHPField(size_t Nr=0, size_t Nz=0,
             Scalar init = Scalar(teuk::zeroC, 0, 0),
             int p = 0, int q = 0)
            : p_(p), q_(q), Nz_(Nz), Nr_(Nr)
            {
                zero_pq_ = Scalar(teuk::zeroC,0,0);
                // resize outer vector
                values_.resize(Nr_);

                // initialize each r–row with Nz Scalars
                for (size_t r = 0; r < Nr_; ++r) {
                    values_[r].resize(Nz_);
                    for (size_t z = 0; z < Nz_; ++z) {
                        values_[r][z] = Scalar(init.value(), p_, q_);
                    }
                }
            }

    /**
      * @brief Fill the field with values computed from a callable f(r,z).
      * @param func  Function taking (r_index, z_index) → Complex
      * The returned Complex is wrapped into a GHPScalar with fixed (p,q).
      */
    void fill(std::function<Complex(size_t, size_t)> func) {
        for (size_t r = 0; r < Nr(); ++r)
            for (size_t z = 0; z < Nz(); ++z)
                values_[r][z] = Scalar(func(r, z), p_, q_);
    }

    const RZBlockScalar& values() const { return values_; }
    void set_values(const RZBlockScalar& vals) { values_ = vals; }

    Scalar& operator()(size_t ir, size_t iz) { return values_[ir][iz]; }
    const Scalar& operator()(size_t ir, size_t iz) const { return values_[ir][iz]; }
    void set(size_t ir, size_t jz, Scalar val) { values_[ir][jz] = val; }

    int p() const { return p_; }
    int q() const { return q_; }
    size_t Nz() const { return Nz_; }
    size_t Nr() const { return Nr_; }

    /**
     * @brief GHP conjugation: swaps weights (p,q) → (q,p) and conjugates components.
     * @return A new GHPField with weights (q,p).
     */
    [[nodiscard]] GHPField conj() const {
        GHPField result(Nz_, Nr_, zero_pq_, q_, p_);
        for(size_t i = 0; i < Nz_; ++i)
            for(size_t j = 0; j < Nr_; ++j)
                result.values_[i][j] = values_[i][j].conj();
        return result;
    }

    /**
   * @brief Apply a GHP spin-boost transformation.
   * @param lambda Complex spin-boost parameter
   * @return New transformed GHPField
   */
    GHPField transform(const Complex& lambda) const {
        GHPField result(Nz_, Nr_,zero_pq_, p_, q_);
        for(size_t i = 0; i < Nz_; ++i)
            for(size_t j = 0; j < Nr_; ++j)
                result.values_[i][j] = math::PowInt(lambda, p_) *
                        math::PowInt(std::conj(lambda), q_) * values_[i][j];
        return result;
    }
//
// Operator overloads
//
    // Element-wise multiplication
     GHPField operator*(const GHPField& other) const {
        if(Nz_ != other.Nz_ || Nr_ != other.Nr_)
            throw std::runtime_error("GHPField dimensions mismatch in multiplication");
        GHPField result(Nz_, Nr_, zero_pq_, p_ + other.p_, q_ + other.q_);
        for(size_t i = 0; i < Nz_; ++i)
            for(size_t j = 0; j < Nr_; ++j)
                result.values_[i][j] = values_[i][j] * other.values_[i][j];
        return result;
    }
    // Element-wise addition (weights must match)
    [[nodiscard]] GHPField operator+(const GHPField& other) const {
        if(Nz_ != other.Nz_ || Nr_ != other.Nr_)
            throw std::runtime_error("GHPField dimensions mismatch in addition");
        if(p_ != other.p_ || q_ != other.q_)
            throw std::runtime_error("GHPField weights mismatch in addition");
        GHPField result(Nz_, Nr_,zero_pq_, p_, q_);
        for(size_t i=0; i<Nz_; ++i)
            for(size_t j=0; j<Nr_; ++j)
                result.values_[i][j] = values_[i][j] + other.values_[i][j];
        return result;
    }

    // Convenience string representation for debugging
    std::string str(size_t ir = -1, size_t jz = -1) const {
        std::ostringstream oss;
        oss << "GHPField(p,q)=(" << p_ << "," << q_ << "), size=(" << Nz_ << "," << Nr_ << ")\n";
        if(ir >= 0 && jz >= 0) {
            const auto& val = values_[ir][jz].value();
            oss << "values_[" << ir << "][" << jz << "] = "
                << std::setprecision(6) << val << "\n";
        }
        return oss.str();
    }

    // Access the underlying 2D array (read-only)
    const RZBlockScalar& data() const { return values_; }

    // Optionally: fill field with a function of (z,r)
    void fill(std::function<Complex(double r, double z)> func,
              const std::vector<double>& r_nodes,
              const std::vector<double>& z_nodes)
    {
        for(size_t i = 0; i < Nr_; ++i) {
            for(size_t j = 0; j < Nz_; ++j) {
                values_[i][j].value() = func(r_nodes[i], z_nodes[j]);
            }
        }
    }
};


enum class GHPCoefficientType {
    kappa, kappap, sigma, sigmap,
    tau, taup, rho, rhop //GHP Covariant coeffs
};

// GHP spin coefficient container
class GHPCoefficients {
private:
    std::unordered_map<GHPCoefficientType, GHPScalar<Complex>> coeffs;  // creates key value pairs

public:
    GHPCoefficients() = default;

    void set(GHPCoefficientType type,  GHPScalar<Complex> value);
    [[nodiscard]]  GHPScalar<Complex> get(GHPCoefficientType type) const;

    // Optional: print all coefficients
    std::string toString() const;
};



struct SpinCoefficientsGHP {
    GHPScalar<Complex> kappa, kappap, sigma, sigmap, tau, taup, rho, rhop; //GHP Covariant coeffs
    GHPScalar<Complex> kappa_bar, kappap_bar, sigma_bar, sigmap_bar, tau_bar, taup_bar, rho_bar, rhop_bar; //GHP Covariant coeffs
    Complex beta, betap, epsilon, epsilonp;
    SpinCoefficientsGHP() = default;

    explicit SpinCoefficientsGHP(const SpinCoefficients& sc_np);

};

#endif // GHZ_NUMERIC_GHPSCALARS_HPP