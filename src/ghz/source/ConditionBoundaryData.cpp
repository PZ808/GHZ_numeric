//
// Created by Peter Zimmerman on 24.03.26.
//
#include "ghz/source/ConditionBoundaryData.hpp"

#include "ghz/source/ConditionBoundaryData.hpp"

#include <sstream>

namespace ghz::source {

    ConditionBoundaryData::ConditionBoundaryData(const BinaryBoundaryDataArchive& archive,
                                                 int m,
                                                 PatchSide patch)
            : meta_(archive.metadata()),
              mode_data_(archive.load_mode_patch(m, patch)),
              m_(m),
              patch_(patch),
              y_nodes_(mode_data_.Y.values)
    {
        if (y_nodes_.size() < 2) {
            throw std::runtime_error(
                    "ConditionBoundaryData: loaded patch has fewer than 2 Y nodes");
        }

        for (std::size_t c = 0; c < boundary_component_count; ++c) {
            y_values_[c] = mode_data_.components[c];
            y_slopes_[c] = pchip_slopes_complex_(y_nodes_, y_values_[c]);
            conditioned_values_[c] = teuk::zeroC;
        }
    }

    void ConditionBoundaryData::set_z_slice(Real z_fixed)
    {
        z_fixed_ = z_fixed;
        has_z_slice_ = true;

        const Real Yq = canonical_y_(z_fixed_);

        const Real ymin = y_nodes_.front();
        const Real ymax = y_nodes_.back();

        if (Yq < ymin || Yq > ymax) {
            std::ostringstream os;
            os << "ConditionBoundaryData::set_z_slice(): requested z=" << z_fixed
               << " maps to Y=" << Yq << " outside stored domain ["
               << ymin << ", " << ymax << "]";
            throw std::runtime_error(os.str());
        }

        for (std::size_t c = 0; c < boundary_component_count; ++c) {
            conditioned_values_[c] = cubic_interp_(y_values_[c], y_slopes_[c], Yq);
        }
    }

    ConditionBoundaryData::Real ConditionBoundaryData::canonical_y_(Real z) const
    {
        const Real Yq = -meta_.rp * z;

        switch (meta_.sym.y_storage) {
            case YStorage::Full:
                return Yq;
            case YStorage::HalfWithEvenReflection:
                return abs_(Yq);
            default:
                throw std::runtime_error(
                        "ConditionBoundaryData::canonical_y_(): unsupported YStorage");
        }
    }

    ConditionBoundaryData::Complex
    ConditionBoundaryData::eval(BoundaryComponent c) const
    {
        if (!has_z_slice_) {
            throw std::runtime_error(
                    "ConditionBoundaryData::eval(): z slice not set; call set_z_slice(z) first");
        }

        return conditioned_values_[static_cast<std::size_t>(c)];
    }

    size_t ConditionBoundaryData::lower_cell_(const std::vector<Real>& grid, Real x)
    {
        if (grid.size() < 2) {
            throw std::runtime_error(
                    "ConditionBoundaryData::lower_cell_(): grid too short");
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

    ConditionBoundaryData::Complex
    ConditionBoundaryData::cubic_interp_(const std::vector<Complex>& f,
                                         const std::vector<Complex>& d,
                                         Real yq) const
    {
        const size_t i0 = lower_cell_(y_nodes_, yq);
        const size_t i1 = i0 + 1;

        const Real y0 = y_nodes_[i0];
        const Real y1 = y_nodes_[i1];
        const Real h  = y1 - y0;

        if (h == Real(0)) {
            return f[i0];
        }

        const Real t   = (yq - y0) / h;
        const Real t2  = t * t;
        const Real t3  = t2 * t;

        const Real h00 =  Real(2) * t3 - Real(3) * t2 + Real(1);
        const Real h10 =         t3 - Real(2) * t2 + t;
        const Real h01 = -Real(2) * t3 + Real(3) * t2;
        const Real h11 =         t3 -         t2;

        return h00 * f[i0] + h10 * h * d[i0]
               + h01 * f[i1] + h11 * h * d[i1];
    }

    bool ConditionBoundaryData::same_sign_nonzero_(Real a, Real b) noexcept
    {
        return ((a > Real(0) && b > Real(0)) ||
                (a < Real(0) && b < Real(0)));
    }

    ConditionBoundaryData::Real ConditionBoundaryData::abs_(Real x) noexcept
    {
        return teuk::Fabs(x);
    }

    std::vector<ConditionBoundaryData::Real>
    ConditionBoundaryData::pchip_slopes_real_(const std::vector<Real>& x,
                                              const std::vector<Real>& y)
    {
        const size_t n = x.size();

        if (n != y.size()) {
            throw std::runtime_error(
                    "ConditionBoundaryData::pchip_slopes_real_(): x/y size mismatch");
        }
        if (n < 2) {
            throw std::runtime_error(
                    "ConditionBoundaryData::pchip_slopes_real_(): need at least 2 nodes");
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
                        "ConditionBoundaryData::pchip_slopes_real_(): x must be strictly increasing");
            }
            delta[i] = (y[i + 1] - y[i]) / h[i];
        }

        for (size_t i = 1; i + 1 < n; ++i) {
            if (same_sign_nonzero_(delta[i - 1], delta[i])) {
                const Real w1 = Real(2) * h[i] + h[i - 1];
                const Real w2 = h[i] + Real(2) * h[i - 1];
                d[i] = (w1 + w2) / (w1 / delta[i - 1] + w2 / delta[i]);
            } else {
                d[i] = Real(0);
            }
        }

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

    std::vector<ConditionBoundaryData::Complex>
    ConditionBoundaryData::pchip_slopes_complex_(const std::vector<Real>& x,
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

} // namespace ghz::source:w
