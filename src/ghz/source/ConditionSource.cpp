#include "ghz/source/ConditionSource.hpp"

#include <sstream>

namespace ghz::source {

    ConditionSource::ConditionSource(const EffectiveSourceArchive& archive,
                                     int m,
                                     PatchSide patch)
            : meta_(archive.metadata()),
              sampler_(meta_, archive.load_mode_patch(m, patch)),
              m_(m),
              patch_(patch),
              x_nodes_(sampler_.mode_data().X.values)
    {
        if (x_nodes_.size() < 2) {
            throw std::runtime_error(
                    "ConditionSource: loaded patch has fewer than 2 X nodes"
            );
        }

        for (auto& v : slice_values_) {
            v.resize(x_nodes_.size(), teuk::zeroC);
        }
        for (auto& d : slice_slopes_) {
            d.resize(x_nodes_.size(), teuk::zeroC);
        }
    }

    void ConditionSource::set_z_slice(Real z_fixed)
    {
        z_fixed_ = z_fixed;
        has_z_slice_ = true;

        const Real Yq = -meta_.rp * z_fixed_;

        // Sample all components on the native X nodes of this patch
        for (size_t c = 0; c < component_count; ++c) {
            const Component comp = static_cast<Component>(c);

            for (size_t ix = 0; ix < x_nodes_.size(); ++ix) {
                slice_values_[c][ix] =
                        sampler_.sample_native(comp, x_nodes_[ix], Yq);
            }

            slice_slopes_[c] = pchip_slopes_complex_(x_nodes_, slice_values_[c]);
        }
    }

    ConditionSource::Real ConditionSource::r_min() const noexcept
    {
        const Real r0 = meta_.rp + meta_.B * x_nodes_.front();
        const Real r1 = meta_.rp + meta_.B * x_nodes_.back();
        return std::min(r0, r1);
    }

    ConditionSource::Real ConditionSource::r_max() const noexcept
    {
        const Real r0 = meta_.rp + meta_.B * x_nodes_.front();
        const Real r1 = meta_.rp + meta_.B * x_nodes_.back();
        return std::max(r0, r1);
    }

    bool ConditionSource::contains_r(Real r) const noexcept
    {
        return (r >= r_min() && r <= r_max());
    }

    ConditionSource::Complex ConditionSource::eval(Component c, Real r) const
    {
        if (!has_z_slice_) {
            throw std::runtime_error(
                    "ConditionSource::eval(): z slice not set; call set_z_slice(z) first"
            );
        }

        const Real Xq = (r - meta_.rp) / meta_.B;

        const Real xmin = x_nodes_.front();
        const Real xmax = x_nodes_.back();

        if (Xq < xmin || Xq > xmax) {
            std::ostringstream os;
            os << "ConditionSource::eval(): requested r=" << r
               << " (X=" << Xq << ") outside patch domain ["
               << xmin << ", " << xmax << "]";
            throw std::runtime_error(os.str());
        }

        const size_t idx = component_index(c);
        return cubic_interp_(slice_values_[idx], slice_slopes_[idx], Xq);
    }

    std::vector<ConditionSource::Complex>
    ConditionSource::eval_on_r_grid(Component c,
                                    const std::vector<Real>& r_nodes) const
    {
        std::vector<Complex> out(r_nodes.size());

#pragma omp parallel for default(none) shared(out, r_nodes, c)
        for (size_t i = 0; i < r_nodes.size(); ++i) {
            out[i] = eval(c, r_nodes[i]);
        }

        return out;
    }

    size_t ConditionSource::lower_cell_(const std::vector<Real>& grid, Real x)
    {
        if (grid.size() < 2) {
            throw std::runtime_error(
                    "ConditionSource::lower_cell_(): grid too short"
            );
        }

        auto it = std::upper_bound(grid.begin(), grid.end(), x);

        if (it == grid.begin()) {
            return 0;
        }

        if (it == grid.end()) {
            return grid.size() - 2;
        }

        return static_cast<size_t>(std::distance(grid.begin(), it) - 1);
    }

    ConditionSource::Complex
    ConditionSource::cubic_interp_(const std::vector<Complex>& f,
                                   const std::vector<Complex>& d,
                                   Real Xq) const
    {
        const size_t i0 = lower_cell_(x_nodes_, Xq);
        const size_t i1 = i0 + 1;

        const Real x0 = x_nodes_[i0];
        const Real x1 = x_nodes_[i1];
        const Real h  = x1 - x0;

        if (h == Real(0)) {
            return f[i0];
        }

        const Real t   = (Xq - x0) / h;
        const Real t2  = t * t;
        const Real t3  = t2 * t;

        const Real h00 =  Real(2) * t3 - Real(3) * t2 + Real(1);
        const Real h10 =         t3 - Real(2) * t2 + t;
        const Real h01 = -Real(2) * t3 + Real(3) * t2;
        const Real h11 =         t3 -         t2;

        return h00 * f[i0] + h10 * h * d[i0]
               + h01 * f[i1] + h11 * h * d[i1];
    }

    bool ConditionSource::same_sign_nonzero_(Real a, Real b) noexcept
    {
        return ((a > Real(0) && b > Real(0)) ||
                (a < Real(0) && b < Real(0)));
    }

    ConditionSource::Real ConditionSource::abs_(Real x) noexcept
    {
        return teuk::Fabs(x);
    }

    std::vector<ConditionSource::Real>
    ConditionSource::pchip_slopes_real_(const std::vector<Real>& x,
                                        const std::vector<Real>& y)
    {
        const size_t n = x.size();

        if (n != y.size()) {
            throw std::runtime_error(
                    "ConditionSource::pchip_slopes_real_(): x/y size mismatch"
            );
        }
        if (n < 2) {
            throw std::runtime_error(
                    "ConditionSource::pchip_slopes_real_(): need at least 2 nodes"
            );
        }

        std::vector<Real> d(n, Real(0));

        if (n == 2) {
            const Real h = x[1] - x[0];
            const Real delta = (y[1] - y[0]) / h;
            d[0] = delta;
            d[1] = delta;
            return d;
        }

        std::vector<Real> h(n - 1), delta(n - 1);
        for (size_t i = 0; i + 1 < n; ++i) {
            h[i] = x[i + 1] - x[i];
            if (!(h[i] > Real(0))) {
                throw std::runtime_error(
                        "ConditionSource::pchip_slopes_real_(): x must be strictly increasing"
                );
            }
            delta[i] = (y[i + 1] - y[i]) / h[i];
        }

        // interior slopes: weighted harmonic mean where appropriate
        for (size_t i = 1; i + 1 < n; ++i) {
            if (same_sign_nonzero_(delta[i - 1], delta[i])) {
                const Real w1 = Real(2) * h[i] + h[i - 1];
                const Real w2 = h[i] + Real(2) * h[i - 1];
                d[i] = (w1 + w2) / (w1 / delta[i - 1] + w2 / delta[i]);
            } else {
                d[i] = Real(0);
            }
        }

        // endpoint slope d[0]
        {
            const Real h0 = h[0];
            const Real h1 = h[1];
            const Real del0 = delta[0];
            const Real del1 = delta[1];

            Real d0 = ((Real(2) * h0 + h1) * del0 - h0 * del1) / (h0 + h1);

            if (!same_sign_nonzero_(d0, del0)) {
                d0 = Real(0);
            } else if (!same_sign_nonzero_(del0, del1) && abs_(d0) > Real(3) * abs_(del0)) {
                d0 = Real(3) * del0;
            }
            d[0] = d0;
        }

        // endpoint slope d[n-1]
        {
            const Real hn2 = h[n - 2];
            const Real hn3 = h[n - 3];
            const Real deln2 = delta[n - 2];
            const Real deln3 = delta[n - 3];

            Real dn = ((Real(2) * hn2 + hn3) * deln2 - hn2 * deln3) / (hn2 + hn3);

            if (!same_sign_nonzero_(dn, deln2)) {
                dn = Real(0);
            } else if (!same_sign_nonzero_(deln2, deln3) && abs_(dn) > Real(3) * abs_(deln2)) {
                dn = Real(3) * deln2;
            }
            d[n - 1] = dn;
        }

        return d;
    }

    std::vector<ConditionSource::Complex>
    ConditionSource::pchip_slopes_complex_(const std::vector<Real>& x,
                                           const std::vector<Complex>& y)
    {
        std::vector<Real> yr(y.size()), yi(y.size());
        for (size_t i = 0; i < y.size(); ++i) {
            yr[i] = y[i].real();
            yi[i] = y[i].imag();
        }

        const std::vector<Real> dr = pchip_slopes_real_(x, yr);
        const std::vector<Real> di = pchip_slopes_real_(x, yi);

        std::vector<Complex> out(y.size());
        for (size_t i = 0; i < y.size(); ++i) {
            out[i] = Complex(dr[i], di[i]);
        }
        return out;
    }

} // namespace ghz::source