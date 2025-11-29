//
// Created by Peter Zimmerman on 28.11.25.
//

#ifndef GHZ_NUMERIC_SPLINES_HPP
#define GHZ_NUMERIC_SPLINES_HPP

#pragma once

#include "GhzTypes.hpp"

#include <vector>
#include <cmath>
#include <cassert>
#include <algorithm>
namespace ghz {

    using Real = teuk::Real;

    class PeriodicSpline {
    public:
        PeriodicSpline() = default;

        // Construct from knots
        PeriodicSpline(const std::vector<Real>& x_input,
                       const std::vector<Real>& y_input) {
            set_points(x_input, y_input);
        }

        [[nodiscard]] Real eval(Real x) const {
            assert(N_ >= 2);

            // wrap x into [x0, x0 + period)
            x = wrap(x);

            // find segment index
            auto it = std::upper_bound(x_.begin(), x_.end(), x);
            std::size_t i = std::clamp(std::size_t(std::distance(x_.begin(), it) - 1), std::size_t(0), N_ - 2);

            const Real h = x_[i + 1] - x_[i];
            const Real a = (x_[i + 1] - x) / h;
            const Real b = (x - x_[i]) / h;

            return a * y_[i] + b * y_[i + 1] +
                   ((a*a*a - a) * M_[i] + (b*b*b - b) * M_[i + 1]) * (h*h) / 6.0;
        }

        Real operator()(Real x) const { return eval(x); }

        [[nodiscard]] Real period() const noexcept { return period_; }
        [[nodiscard]] std::size_t size() const noexcept { return N_; }
        [[nodiscard]] const std::vector<Real>& x_grid() const noexcept { return x_; }
        [[nodiscard]] const std::vector<Real>& y_grid() const noexcept { return y_; }

    private:
        std::size_t N_ = 0;
        Real period_ = 0.0;

        std::vector<Real> x_, y_, M_;

        [[nodiscard]] Real wrap(Real x) const noexcept {
            const Real x0 = x_.front();
            Real xp = fmod(x - x0, period_);
            if (xp < 0) xp += period_;
            return x0 + xp;
        }

        void set_points(const std::vector<Real>& xin,
                        const std::vector<Real>& yin) {
            assert(xin.size() == yin.size() && xin.size() >= 2);

            N_ = xin.size();
            x_ = xin;
            y_ = yin;
            period_ = x_.back() - x_.front();

            M_.resize(N_, 0.0);
            if (N_ == 2) return; // linear case

            // Compute h and system coefficients
            std::vector<Real> h(N_ - 1);
            for (std::size_t i = 0; i < N_ - 1; ++i) h[i] = x_[i + 1] - x_[i];

            std::vector<Real> a(N_), b(N_), c(N_), d(N_);
            for (std::size_t i = 1; i < N_ - 1; ++i) {
                a[i] = h[i - 1];
                b[i] = 2 * (h[i - 1] + h[i]);
                c[i] = h[i];
                d[i] = 6 * ((y_[i + 1] - y_[i]) / h[i] - (y_[i] - y_[i - 1]) / h[i - 1]);
            }

            // Solve reduced tridiagonal system for periodic spline (Sherman-Morrison)
            const std::size_t n = N_ - 2;
            std::vector<Real> cp(n), dp(n);

            cp[0] = c[1] / b[1];
            dp[0] = d[1] / b[1];
            for (std::size_t i = 1; i < n; ++i) {
                const Real denom = b[i + 1] - a[i + 1] * cp[i - 1];
                cp[i] = c[i + 1] / denom;
                dp[i] = (d[i + 1] - a[i + 1] * dp[i - 1]) / denom;
            }

            std::vector<Real> M_small(n);
            M_small[n - 1] = dp[n - 1];
            for (int i = int(n) - 2; i >= 0; --i)
                M_small[i] = dp[i] - cp[i] * M_small[i + 1];

            for (std::size_t i = 1; i < N_ - 1; ++i)
                M_[i] = M_small[i - 1];

            M_[0] = M_[N_ - 1] = M_small[0]; // periodic
        }
    };

} // namespace ghz
#endif //GHZ_NUMERIC_SPLINES_HPP
