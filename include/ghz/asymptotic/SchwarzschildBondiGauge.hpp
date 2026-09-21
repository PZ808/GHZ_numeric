#ifndef GHZ_ASYMPTOTIC_SCHWARZSCHILD_BONDI_GAUGE_HPP
#define GHZ_ASYMPTOTIC_SCHWARZSCHILD_BONDI_GAUGE_HPP

#include "ghz/asymptotic/SchwarzschildPsi0Modes.hpp"
#include "ghz/orbit/KerrOrbit.hpp"

#include <vector>

namespace ghz::asymptotic {

struct HeldGaugeVectorMode {
    teuk::Complex l;
    teuk::Complex n;
    teuk::Complex m;
};

struct HeldGaugeVectorGrid {
    std::vector<teuk::Complex> l;
    std::vector<teuk::Complex> n;
    std::vector<teuk::Complex> m;
};

struct NPGaugeVectorMode {
    teuk::Complex l;
    teuk::Complex n;
    teuk::Complex m;
    teuk::Complex mbar;
};

/** Schwarzschild non-stationary Bondi gauge vector reconstructed from psi0^5o.
 *
 * The input coefficients and all returned amplitudes use the mostly-minus
 * signature.  Harmonic phase factors exp(-i omega u + i m phi) are omitted
 * from coefficient-returning methods.
 */
class SchwarzschildBondiGauge {
public:
    explicit SchwarzschildBondiGauge(const SchwarzschildPsi0Modes& psi0)
        : psi0_(psi0) {}

    [[nodiscard]] HeldGaugeVectorMode held_mode(
        int ell, int m, teuk::Real omega) const;

    /** Sum ell at fixed m after removing each component's regular pole factor.
     * l and n have spin weight 0; m has spin weight +1.
     */
    [[nodiscard]] HeldGaugeVectorGrid summed_reduced_held_mode(
        int m, teuk::Real omega, const std::vector<teuk::Real>& z) const;

    [[nodiscard]] NPGaugeVectorMode np_mode_at_radius(
        int ell, int m, teuk::Real omega, teuk::Real radius) const;

    [[nodiscard]] teuk::Complex u_dot_zeta_amplitude(
        int ell, int m, teuk::Real omega,
        const orbit::KerrCircularEquatorialOrbit& orbit) const;

    /** Pullback u^a u^b (L_zeta g)_ab for a circular geodesic mode.
     *
     * For a geodesic this is 2 u^a nabla_a(u.zeta).  The modal derivative is
     * i(m u^phi - omega u^u), which vanishes for omega=m Omega_phi.
     */
    [[nodiscard]] teuk::Complex lie_derivative_uu(
        int ell, int m, teuk::Real omega,
        const orbit::KerrCircularEquatorialOrbit& orbit,
        teuk::Real outgoing_u = 0, teuk::Real phi = 0) const;

private:
    const SchwarzschildPsi0Modes& psi0_;
};

} // namespace ghz::asymptotic

#endif
