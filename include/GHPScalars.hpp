//
// Created by Peter Zimmerman on 25.10.25.
//

#ifndef GHZ_NUMERIC_GHPSCALARS_HPP
#define GHZ_NUMERIC_GHPSCALARS_HPP
#include <complex>
#include "SpinCoeffsNP.hpp"

//using Complex = std::complex<double>;

class GHPScalar {
public:
    using Complex = teuk::Complex;

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
    GHPScalar kappa, kappap, sigma, sigmap, tau, taup, rho, rhop;
    SpinCoefficientsGHP() = default;

    explicit SpinCoefficientsGHP(const SpinCoefficients& sc_np);

    // Apply the GHP prime operation: swap tetrad roles
    // [[nodiscard]] SpinCoefficientsGHP prime() const;
};



#endif //GHZ_NUMERIC_GHPSCALARS_HPP
