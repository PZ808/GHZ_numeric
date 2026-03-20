//
// Created by Peter Zimmerman on 28.11.25.
//

#ifndef GHZ_NUMERIC_SPLINES_HPP
#define GHZ_NUMERIC_SPLINES_HPP

#pragma once

#include "ghz/core/GhzTypes.hpp"

#include <vector>
#include <cmath>
#include <cassert>
#include <algorithm>

namespace orbit {

    using Real = teuk::Real;

    class PeriodicSpline {
    public:
        PeriodicSpline() = default;

        // Construct from knots
        PeriodicSpline(const std::vector<Real>& x_input,
                       const std::vector<Real>& y_input) {
            set_points(x_input, y_input);
        }
        PeriodicSpline(const std::vector<Real>& x_input,
                       const std::vector<Real>& y_input,
                       Real period_input) {
            set_points(x_input, y_input, period_input);
        }

        Real eval(Real x) const {
            assert(N_ >= 2);

            x = wrap(x);

            // periodic last interval: [x_[N_-1], x_[0] + period_)
            if (x >= x_.back()) {
                const Real xL = x_.back();
                const Real xR = x_.front() + period_;
                const Real h  = xR - xL;
                const Real a  = (xR - x) / h;
                const Real b  = (x - xL) / h;

                return a * y_.back() + b * y_.front()
                       + ((a*a*a - a) * M_.back() + (b*b*b - b) * M_.front()) * (h*h) / 6.0;
            }

            auto it = std::upper_bound(x_.begin(), x_.end(), x);
            std::size_t i = std::size_t(std::distance(x_.begin(), it) - 1);

            const Real h = x_[i + 1] - x_[i];
            const Real a = (x_[i + 1] - x) / h;
            const Real b = (x - x_[i]) / h;

            return a * y_[i] + b * y_[i + 1]
                   + ((a*a*a - a) * M_[i] + (b*b*b - b) * M_[i + 1]) * (h*h) / 6.0;
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

            // Strictly increasing knots, no duplicated endpoint.
            for (std::size_t i = 1; i < N_; ++i) {
                assert(x_[i] > x_[i - 1]);
            }

            // For a periodic grid without duplicating the endpoint,
            // the period is one extra spacing beyond x_.back().
            period_ = (x_.back() - x_.front()) + (x_[1] - x_[0]);

            M_.assign(N_, Real(0));

            if (N_ == 2) {
                // Degenerate linear periodic-ish fallback.
                return;
            }

            // Interval lengths h_i for intervals [x_i, x_{i+1}],
            // with the last one being [x_{N-1}, x_0 + period_].
            std::vector<Real> h(N_);
            for (std::size_t i = 0; i < N_ - 1; ++i) {
                h[i] = x_[i + 1] - x_[i];
            }
            h[N_ - 1] = (x_.front() + period_) - x_.back();

            // Build the full cyclic system A M = rhs
            std::vector<std::vector<Real>> A(N_, std::vector<Real>(N_, Real(0)));
            std::vector<Real> rhs(N_, Real(0));

            auto mod = [this](long long i) -> std::size_t {
                const long long n = static_cast<long long>(N_);
                long long r = i % n;
                if (r < 0) r += n;
                return static_cast<std::size_t>(r);
            };

            for (std::size_t i = 0; i < N_; ++i) {
                const std::size_t im1 = mod(static_cast<long long>(i) - 1);
                const std::size_t ip1 = mod(static_cast<long long>(i) + 1);

                A[i][im1] = h[im1];
                A[i][i]   = Real(2) * (h[im1] + h[i]);
                A[i][ip1] = h[i];

                const Real slope_right = (y_[ip1] - y_[i])   / h[i];
                const Real slope_left  = (y_[i]   - y_[im1]) / h[im1];
                rhs[i] = Real(6) * (slope_right - slope_left);
            }

            // Solve dense linear system with Gaussian elimination + partial pivoting.
            for (std::size_t k = 0; k < N_; ++k) {
                // Pivot
                std::size_t piv = k;
                Real max_abs = std::abs(A[k][k]);
                for (std::size_t i = k + 1; i < N_; ++i) {
                    Real v = std::abs(A[i][k]);
                    if (v > max_abs) {
                        max_abs = v;
                        piv = i;
                    }
                }
                assert(max_abs > Real(0) && "PeriodicSpline: singular linear system");

                if (piv != k) {
                    std::swap(A[k], A[piv]);
                    std::swap(rhs[k], rhs[piv]);
                }

                // Eliminate
                for (std::size_t i = k + 1; i < N_; ++i) {
                    const Real factor = A[i][k] / A[k][k];
                    if (factor == Real(0)) continue;

                    A[i][k] = Real(0);
                    for (std::size_t j = k + 1; j < N_; ++j) {
                        A[i][j] -= factor * A[k][j];
                    }
                    rhs[i] -= factor * rhs[k];
                }
            }

            // Back substitution
            for (long long i = static_cast<long long>(N_) - 1; i >= 0; --i) {
                Real sum = rhs[static_cast<std::size_t>(i)];
                for (std::size_t j = static_cast<std::size_t>(i) + 1; j < N_; ++j) {
                    sum -= A[static_cast<std::size_t>(i)][j] * M_[j];
                }
                M_[static_cast<std::size_t>(i)] = sum / A[static_cast<std::size_t>(i)][static_cast<std::size_t>(i)];
            }
        }

        void set_points(const std::vector<Real>& xin,
                        const std::vector<Real>& yin,
                        Real period_input) {
            assert(xin.size() == yin.size() && xin.size() >= 2);

            N_ = xin.size();
            x_ = xin;
            y_ = yin;
            period_ = period_input;

            M_.assign(N_, 0.0);
            if (N_ == 2) return;

            // ... same cyclic-system code as before ...
            //
            // Interval lengths h_i for intervals [x_i, x_{i+1}],
            // with the last one being [x_{N-1}, x_0 + period_].
            std::vector<Real> h(N_);
            for (std::size_t i = 0; i < N_ - 1; ++i) {
                h[i] = x_[i + 1] - x_[i];
            }
            h[N_ - 1] = (x_.front() + period_) - x_.back();

            // Build the full cyclic system A M = rhs
            std::vector<std::vector<Real>> A(N_, std::vector<Real>(N_, Real(0)));
            std::vector<Real> rhs(N_, Real(0));

            auto mod = [this](long long i) -> std::size_t {
                const long long n = static_cast<long long>(N_);
                long long r = i % n;
                if (r < 0) r += n;
                return static_cast<std::size_t>(r);
            };

            for (std::size_t i = 0; i < N_; ++i) {
                const std::size_t im1 = mod(static_cast<long long>(i) - 1);
                const std::size_t ip1 = mod(static_cast<long long>(i) + 1);

                A[i][im1] = h[im1];
                A[i][i]   = Real(2) * (h[im1] + h[i]);
                A[i][ip1] = h[i];

                const Real slope_right = (y_[ip1] - y_[i])   / h[i];
                const Real slope_left  = (y_[i]   - y_[im1]) / h[im1];
                rhs[i] = Real(6) * (slope_right - slope_left);
            }

            // Solve dense linear system with Gaussian elimination + partial pivoting.
            for (std::size_t k = 0; k < N_; ++k) {
                // Pivot
                std::size_t piv = k;
                Real max_abs = std::abs(A[k][k]);
                for (std::size_t i = k + 1; i < N_; ++i) {
                    Real v = std::abs(A[i][k]);
                    if (v > max_abs) {
                        max_abs = v;
                        piv = i;
                    }
                }
                assert(max_abs > Real(0) && "PeriodicSpline: singular linear system");

                if (piv != k) {
                    std::swap(A[k], A[piv]);
                    std::swap(rhs[k], rhs[piv]);
                }

                // Eliminate
                for (std::size_t i = k + 1; i < N_; ++i) {
                    const Real factor = A[i][k] / A[k][k];
                    if (factor == Real(0)) continue;

                    A[i][k] = Real(0);
                    for (std::size_t j = k + 1; j < N_; ++j) {
                        A[i][j] -= factor * A[k][j];
                    }
                    rhs[i] -= factor * rhs[k];
                }
            }

            // Back substitution
            for (long long i = static_cast<long long>(N_) - 1; i >= 0; --i) {
                Real sum = rhs[static_cast<std::size_t>(i)];
                for (std::size_t j = static_cast<std::size_t>(i) + 1; j < N_; ++j) {
                    sum -= A[static_cast<std::size_t>(i)][j] * M_[j];
                }
                M_[static_cast<std::size_t>(i)] = sum / A[static_cast<std::size_t>(i)][static_cast<std::size_t>(i)];
            }


        }
    };

} // namespace orbit
#endif //GHZ_NUMERIC_SPLINES_HPP
