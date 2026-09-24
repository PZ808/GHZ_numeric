//
// TeukolskySource.cpp
//

#include "ghz/source/TeukolskySource.hpp"

#include "ghz/ghp/GHPScalars.hpp"

#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

namespace ghz::source {
namespace {

    using Complex = teuk::Complex;
    using Real = teuk::Real;
    using Scalar = ghp::GHPScalar<Complex>;

    void require_same_grid(const TeukolskyField& a,
                           const TeukolskyField& b,
                           const char* where) {
        if (a.Nr() != b.Nr() || a.Nz() != b.Nz()) {
            std::ostringstream os;
            os << where << ": field grid mismatch";
            throw std::runtime_error(os.str());
        }
    }

    void require_same_modes(const TeukolskyField& a,
                            const TeukolskyField& b,
                            const char* where) {
        if (a.m() != b.m() || a.kr() != b.kr() || a.kz() != b.kz()) {
            std::ostringstream os;
            os << where << ": source mode mismatch";
            throw std::runtime_error(os.str());
        }
        if (a.has_omega_mk() != b.has_omega_mk()) {
            std::ostringstream os;
            os << where << ": omega metadata mismatch";
            throw std::runtime_error(os.str());
        }
        if (a.has_omega_mk() && a.omega_mk() != b.omega_mk()) {
            std::ostringstream os;
            os << where << ": omega value mismatch";
            throw std::runtime_error(os.str());
        }
    }

    void require_source_compatible(const TeukolskyField& ref,
                                   const TeukolskyField& f,
                                   const char* where) {
        require_same_grid(ref, f, where);
        require_same_modes(ref, f, where);
    }

    void require_coeff_compatible(const TeukolskyField& ref,
                                  const TeukolskyField& coeff,
                                  const char* where) {
        require_same_grid(ref, coeff, where);
    }

    TeukolskyField make_like(const TeukolskyField& ref,
                             int p,
                             int q) {
        TeukolskyField out(ref.Nr(),
                           ref.Nz(),
                           ref.modes(),
                           Scalar(teuk::zeroC, p, q),
                           p,
                           q);
        if (ref.has_omega_mk()) {
            out.set_omega_mk(ref.omega_mk());
        }
        return out;
    }

    void require_weight(const TeukolskyField& field,
                        int p,
                        int q,
                        const char* where) {
        if (field.p() != p || field.q() != q) {
            std::ostringstream os;
            os << where << ": field has weights (" << field.p() << ","
               << field.q() << "), expected (" << p << "," << q << ")";
            throw std::runtime_error(os.str());
        }
    }

    TeukolskyField apply_thorn(const TeukolskyField& in,
                               const TeukolskySourceOperators& ops) {
        if (!ops.thorn) {
            throw std::runtime_error("apply_s0_teukolsky_source: missing thorn operator");
        }
        TeukolskyField out = make_like(in, in.p() + 1, in.q() + 1);
        ops.thorn(in, out);
        return out;
    }

    TeukolskyField apply_eth(const TeukolskyField& in,
                             const TeukolskySourceOperators& ops) {
        if (!ops.eth) {
            throw std::runtime_error("apply_s0_teukolsky_source: missing eth operator");
        }
        TeukolskyField out = make_like(in, in.p(), in.q() - 2);
        ops.eth(in, out);
        return out;
    }

    TeukolskyField coeff_times(const TeukolskyField& coeff,
                               const TeukolskyField& field,
                               const char* where) {
        require_coeff_compatible(field, coeff, where);
        TeukolskyField out =
                make_like(field, coeff.p() + field.p(), coeff.q() + field.q());

        for (size_t ir = 0; ir < field.Nr(); ++ir) {
            for (size_t iz = 0; iz < field.Nz(); ++iz) {
                out(ir, iz) = coeff(ir, iz) * field(ir, iz);
            }
        }

        return out;
    }

    TeukolskyField linear_combination(
            std::vector<std::pair<Real, const TeukolskyField*>> terms,
            const char* where) {
        if (terms.empty()) {
            throw std::runtime_error("linear_combination: no terms");
        }

        const auto& ref = *terms.front().second;
        TeukolskyField out = make_like(ref, ref.p(), ref.q());

        for (const auto& [factor, term] : terms) {
            require_source_compatible(ref, *term, where);
            require_weight(*term, ref.p(), ref.q(), where);
        }

        for (size_t ir = 0; ir < ref.Nr(); ++ir) {
            for (size_t iz = 0; iz < ref.Nz(); ++iz) {
                Complex value = teuk::zeroC;
                for (const auto& [factor, term] : terms) {
                    value += factor * (*term)(ir, iz).value();
                }
                out(ir, iz) = Scalar(value, ref.p(), ref.q());
            }
        }

        return out;
    }

    TeukolskyField scale(const TeukolskyField& field,
                         Real factor) {
        TeukolskyField out = make_like(field, field.p(), field.q());
        for (size_t ir = 0; ir < field.Nr(); ++ir) {
            for (size_t iz = 0; iz < field.Nz(); ++iz) {
                out(ir, iz) = Scalar(factor * field(ir, iz).value(),
                                     field.p(),
                                     field.q());
            }
        }
        return out;
    }

} // namespace

TeukolskyField apply_s0_teukolsky_source(
        const TeukolskySourceComponents& components,
        const TeukolskySourceCoefficients& coeffs,
        const TeukolskySourceOperators& ops) {
    const auto& Tll = components.Tll;
    const auto& Tlm = components.Tlm_sym;
    const auto& Tmm = components.Tmm;

    require_source_compatible(Tlm, Tll, "apply_s0_teukolsky_source");
    require_source_compatible(Tlm, Tmm, "apply_s0_teukolsky_source");
    require_coeff_compatible(Tlm, coeffs.rho, "apply_s0_teukolsky_source: rho");
    require_coeff_compatible(Tlm, coeffs.rho_bar, "apply_s0_teukolsky_source: rho_bar");
    require_coeff_compatible(Tlm, coeffs.tau, "apply_s0_teukolsky_source: tau");
    require_coeff_compatible(Tlm,
                             coeffs.tau_prime_bar,
                             "apply_s0_teukolsky_source: tau_prime_bar");

    const TeukolskyField thorn_Tlm = apply_thorn(Tlm, ops);
    const TeukolskyField rhobar_Tlm =
            coeff_times(coeffs.rho_bar, Tlm, "rho_bar T_(lm)");
    const TeukolskyField eth_Tll = apply_eth(Tll, ops);
    const TeukolskyField taupbar_Tll =
            coeff_times(coeffs.tau_prime_bar, Tll, "tau_prime_bar T_ll");

    const TeukolskyField first_bracket = linear_combination(
            {{Real(1), &thorn_Tlm},
             {Real(-2), &rhobar_Tlm},
             {Real(-1), &eth_Tll},
             {Real(1), &taupbar_Tll}},
            "first Teukolsky source bracket");

    const TeukolskyField eth_Tlm = apply_eth(Tlm, ops);
    const TeukolskyField taupbar_Tlm =
            coeff_times(coeffs.tau_prime_bar, Tlm, "tau_prime_bar T_(lm)");
    const TeukolskyField thorn_Tmm = apply_thorn(Tmm, ops);
    const TeukolskyField rhobar_Tmm =
            coeff_times(coeffs.rho_bar, Tmm, "rho_bar T_mm");

    const TeukolskyField second_bracket = linear_combination(
            {{Real(1), &eth_Tlm},
             {Real(-2), &taupbar_Tlm},
             {Real(-1), &thorn_Tmm},
             {Real(1), &rhobar_Tmm}},
            "second Teukolsky source bracket");

    const TeukolskyField eth_first = apply_eth(first_bracket, ops);
    const TeukolskyField taupbar_first =
            coeff_times(coeffs.tau_prime_bar,
                         first_bracket,
                         "tau_prime_bar first bracket");
    const TeukolskyField tau_first =
            coeff_times(coeffs.tau, first_bracket, "tau first bracket");

    const TeukolskyField outer_first = linear_combination(
            {{Real(1), &eth_first},
             {Real(-1), &taupbar_first},
             {Real(-4), &tau_first}},
            "outer first Teukolsky source term");

    const TeukolskyField thorn_second = apply_thorn(second_bracket, ops);
    const TeukolskyField rho_second =
            coeff_times(coeffs.rho, second_bracket, "rho second bracket");
    const TeukolskyField rhobar_second =
            coeff_times(coeffs.rho_bar,
                         second_bracket,
                         "rho_bar second bracket");

    const TeukolskyField outer_second = linear_combination(
            {{Real(1), &thorn_second},
             {Real(-4), &rho_second},
             {Real(-1), &rhobar_second}},
            "outer second Teukolsky source term");

    const TeukolskyField total = linear_combination(
            {{Real(1), &outer_first}, {Real(1), &outer_second}},
            "full Teukolsky source");

    return scale(total, Real(0.5));
}

} // namespace ghz::source
