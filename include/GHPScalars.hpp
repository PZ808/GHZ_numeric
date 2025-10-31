//
// Created by Peter Zimmerman on 25.10.25.
//

#ifndef GHZ_NUMERIC_GHPSCALARS_HPP
#define GHZ_NUMERIC_GHPSCALARS_HPP
#include <complex>
#include <iomanip>
#include "SpinCoeffsNP.hpp"

template <typename ComplexT = teuk::Complex>
class GHPScalar {
public:
    using Complex = ComplexT;

private:
    Complex value_;
    int p_;
    int q_;

public:
    explicit GHPScalar(Complex val = 0.0, int p = 0, int q = 0)
            : value_(val), p_(p), q_(q) {}

    [[nodiscard]] Complex value() const { return value_; }
    [[nodiscard]] int p() const { return p_; }
    [[nodiscard]] int q() const { return q_; }

    void setValue(Complex val) { value_ = val; }

    // conjugation flips weights (p,q) -> (q,p)
    [[nodiscard]] GHPScalar conj() const {
        return GHPScalar(std::conj(value_), q_, p_);
    }
    // Apply the GHP prime operation: swap tetrad roles
    [[nodiscard]] GHPScalar prime() const;

    // multiplication combines weights
    [[nodiscard]] GHPScalar operator*(const GHPScalar& other) const {
        return GHPScalar(value_ * other.value_, p_ + other.p_, q_ + other.q_);
    }

    [[nodiscard]] GHPScalar operator/(const GHPScalar& other) const {
        return GHPScalar(value_ / other.value_, p_ - other.p_, q_ - other.q_);
    }

    [[nodiscard]] GHPScalar operator+(const GHPScalar& other) const {
        // addition only makes sense if weights match
        if (p_ != other.p_ || q_ != other.q_) {
            throw std::runtime_error("Incompatible GHP weights in addition");
        }
        return GHPScalar(value_ + other.value_, p_, q_);
    }

    // spin-boost transformation
    [[nodiscard]] GHPScalar transform(const Complex& lambda) const {
        Complex newVal = std::pow(lambda, p_) * std::pow(std::conj(lambda), q_) * value_;
        return GHPScalar(newVal, p_, q_);
    }

    std::string str() const {
        return "val=" + std::to_string(value_.real()) +
               (value_.imag() >= 0 ? "+" : "") + std::to_string(value_.imag()) + "i"
               + " (p,q)=(" + std::to_string(p_) + "," + std::to_string(q_) + ")";
    }
};

struct SpinCoefficientsGHP {
    GHPScalar<Complex> kappa, kappap, sigma, sigmap, tau, taup, rho, rhop;
    Complex beta, betap, epsilon, epsilonp;
    SpinCoefficientsGHP() = default;

    explicit SpinCoefficientsGHP(const SpinCoefficients& sc_np);

    // Apply the GHP prime operation: swap tetrad roles
    // [[nodiscard]] SpinCoefficientsGHP prime() const;
};

class GHPField {
public:
    using Complex = teuk::Complex;

private:
    std::vector<std::vector<Complex>> values_; // [Nz][Nphi]
    int p_;
    int q_;
    int Nz_;
    int Nphi_;

public:
    // Constructor: create a Nz x Nphi field, optionally with an initial value
    GHPField(int Nz, int Nphi, Complex init = 0.0, int p = 0, int q = 0)
            : values_(Nz, std::vector<Complex>(Nphi, init)), p_(p), q_(q), Nz_(Nz), Nphi_(Nphi) {}


    [[nodiscard]] const std::vector<std::vector<Complex>>& values() const { return values_; }
    void set_values(const std::vector<std::vector<Complex>>& vals) { values_ = vals; }

    [[nodiscard]] Complex operator()(int iz, int jphi) const { return values_[iz][jphi]; }
    void set(int iz, int jphi, Complex val) { values_[iz][jphi] = val; }

    [[nodiscard]] int p() const { return p_; }
    [[nodiscard]] int q() const { return q_; }
    [[nodiscard]] int Nz() const { return Nz_; }
    [[nodiscard]] int Nphi() const { return Nphi_; }

    // Element-wise conjugation (swaps weights)
    [[nodiscard]] GHPField conj() const {
        GHPField result(Nz_, Nphi_, 0.0, q_, p_);
        for(int i = 0; i < Nz_; ++i)
            for(int j = 0; j < Nphi_; ++j)
                result.values_[i][j] = std::conj(values_[i][j]);
        return result;
    }

    // Prime operation (swaps tetrad roles, may optionally transform values)
    [[nodiscard]] GHPField prime() const {
        // Default: just swap weights
        GHPField result(Nz_, Nphi_, 0.0, -p_, -q_);
        for(int i = 0; i < Nz_; ++i)
            for(int j = 0; j < Nphi_; ++j)
                result.values_[i][j] = values_[i][j]; // optionally add more rules here
        return result;
    }

    // Element-wise multiplication
    [[nodiscard]] GHPField operator*(const GHPField& other) const {
        if(Nz_ != other.Nz_ || Nphi_ != other.Nphi_)
            throw std::runtime_error("GHPField dimensions mismatch in multiplication");
        GHPField result(Nz_, Nphi_, 0.0, p_ + other.p_, q_ + other.q_);
        for(int i = 0; i < Nz_; ++i)
            for(int j = 0; j < Nphi_; ++j)
                result.values_[i][j] = values_[i][j] * other.values_[i][j];
        return result;
    }

    // Element-wise addition (weights must match)
    [[nodiscard]] GHPField operator+(const GHPField& other) const {
        if(Nz_ != other.Nz_ || Nphi_ != other.Nphi_)
            throw std::runtime_error("GHPField dimensions mismatch in addition");
        if(p_ != other.p_ || q_ != other.q_)
            throw std::runtime_error("GHPField weights mismatch in addition");
        GHPField result(Nz_, Nphi_, 0.0, p_, q_);
        for(int i = 0; i < Nz_; ++i)
            for(int j = 0; j < Nphi_; ++j)
                result.values_[i][j] = values_[i][j] + other.values_[i][j];
        return result;
    }

    // Spin-boost transformation (element-wise)
    [[nodiscard]] GHPField transform(const Complex& lambda) const {
        GHPField result(Nz_, Nphi_, 0.0, p_, q_);
        for(int i = 0; i < Nz_; ++i)
            for(int j = 0; j < Nphi_; ++j)
                result.values_[i][j] = std::pow(lambda, p_) * std::pow(std::conj(lambda), q_) * values_[i][j];
        return result;
    }

    // Convenience string representation for debugging
    std::string str(int iz = -1, int jphi = -1) const {
        std::ostringstream oss;
        oss << "GHPField(p,q)=(" << p_ << "," << q_ << "), size=(" << Nz_ << "," << Nphi_ << ")\n";
        if(iz >= 0 && jphi >= 0) {
            const auto& val = values_[iz][jphi];
            oss << "values_[" << iz << "][" << jphi << "] = "
                << std::setprecision(6) << val << "\n";
        }
        return oss.str();
    }

    // Access the underlying 2D array (read-only)
    [[nodiscard]] const std::vector<std::vector<Complex>>& data() const { return values_; }

    // Optionally: fill field with a function of (z,phi)
    void fill(std::function<Complex(double z, double phi)> func,
              const std::vector<double>& z_nodes,
              const std::vector<double>& phi_nodes)
    {
        for(int i = 0; i < Nz_; ++i) {
            for(int j = 0; j < Nphi_; ++j) {
                values_[i][j] = func(z_nodes[i], phi_nodes[j]);
            }
        }
    }
};


#endif //GHZ_NUMERIC_GHPSCALARS_HPP
