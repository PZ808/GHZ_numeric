// Diagnostic of the asymptotic homogeneous adjustment, not the full sourced
// metric. Radial derivatives are analytic jets; angular derivatives use the
// production pole-factorized Held m-mode operators. No ell projection is used
// in the solve.
#include "ghz/asymptotic/BondiHeldMetricGrid.hpp"
#include "ghz/asymptotic/BondiEllMBanded.hpp"
#include "ghz/asymptotic/SchwarzschildPsi0Modes.hpp"
#include "ghz/geom/KerrMetric.hpp"
#include "ghz/geom/KinnersleyTetrad.hpp"
#include "ghz/ghp/HeldScalars.hpp"
#include "ghz/spectral/SpectralCoordinateMaps.hpp"
#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
using R = teuk::Real;
using C = teuk::Complex;
using V = std::vector<C>;
const C imag_unit(0, 1);
using Harm = ghz::asymptotic::SchwarzschildPsi0Modes;
// Only derivatives through order two are needed by Sdagger.
struct RadialJet {
  C v = 0, d = 0, dd = 0;
};
RadialJet operator+(RadialJet x, RadialJet y) {
  return {x.v + y.v, x.d + y.d, x.dd + y.dd};
}
RadialJet operator-(RadialJet x, RadialJet y) {
  return {x.v - y.v, x.d - y.d, x.dd - y.dd};
}
RadialJet operator*(RadialJet x, RadialJet y) {
  return {x.v * y.v, x.d * y.v + x.v * y.d,
          x.dd * y.v + C(2) * x.d * y.d + x.v * y.dd};
}
RadialJet operator*(C c, RadialJet x) { return {c * x.v, c * x.d, c * x.dd}; }
struct ReducedField {
  int p, q, m;
  R w;
  std::vector<RadialJet> x;
  int s() const { return (p - q) / 2; }
};
ReducedField operator*(C c, ReducedField x) {
  for (auto &v : x.x)
    v = c * v;
  return x;
}
ReducedField operator+(ReducedField x, const ReducedField &y) {
  if (x.m != y.m || x.s() != y.s())
    throw std::runtime_error("field sum spin/m mismatch");
  for (size_t j = 0; j < x.x.size(); ++j)
    x.x[j] = x.x[j] + y.x[j];
  return x;
}
ReducedField operator-(ReducedField x, const ReducedField &y) {
  return x + C(-1) * y;
}
ReducedField conj(ReducedField x) {
  std::swap(x.p, x.q);
  x.m = -x.m;
  x.w = -x.w;
  for (auto &v : x.x)
    v = {std::conj(v.v), std::conj(v.d), std::conj(v.dd)};
  return x;
}
ReducedField weights(ReducedField x, int p, int q) {
  x.p = p;
  x.q = q;
  return x;
}
class KerrDiagnostic {
public:
  const spectral::SpectralDiffer &diff;
  const KinnersleyHeldOperators<OutgoingCoords, spectral::SpectralDiffer> &op;
  R a, M = 1;
  KerrDiagnostic(const spectral::SpectralDiffer &d,
                 const KinnersleyHeldOperators<OutgoingCoords,
                                               spectral::SpectralDiffer> &o,
                 R aa)
      : diff(d), op(o), a(aa) {}
  ReducedField field(int p, int q, int m, R w, C v = 0) const {
    return {p, q, m, w, std::vector<RadialJet>(diff.Nz(), {v, 0, 0})};
  }
  ReducedField mul(const ReducedField &x, const ReducedField &y) const {
    ReducedField out = field(x.p + y.p, x.q + y.q, x.m + y.m, x.w + y.w);
    int s = x.s() + y.s(), m = x.m + y.m;
    int ep =
        (std::abs(x.m + x.s()) + std::abs(y.m + y.s()) - std::abs(m + s)) / 2;
    int em =
        (std::abs(x.m - x.s()) + std::abs(y.m - y.s()) - std::abs(m - s)) / 2;
    for (size_t j = 0; j < diff.Nz(); ++j) {
      R z = diff.lgl_nodes()[j];
      out.x[j] =
          C(std::pow(1 - z, ep) * std::pow(1 + z, em)) * (x.x[j] * y.x[j]);
    }
    return out;
  }
  ReducedField h(ReducedField x, bool up = true) const {
    ReducedField out = field(x.p - (up ? 0 : 2), x.q - (up ? 2 : 0), x.m, x.w);
    for (int k = 0; k < 3; ++k) {
      spectral::SpectralGHPVectorized in(
          1, diff.Nz(), {x.m, 0, 0}, GHPScalar<C>(C(0), x.p, x.q), x.p, x.q);
      spectral::SpectralGHPVectorized dest(1, diff.Nz(), {x.m, 0, 0},
                                           GHPScalar<C>(C(0), out.p, out.q),
                                           out.p, out.q);
      in.set_omega_mk(x.w);
      dest.set_omega_mk(x.w);
      for (size_t j = 0; j < diff.Nz(); ++j)
        in(0, j).value() = k == 0 ? x.x[j].v : k == 1 ? x.x[j].d : x.x[j].dd;
      if (up)
        op.edthHRed_inplace(in, dest);
      else
        op.edthBarHRed_inplace(in, dest);
      for (size_t j = 0; j < diff.Nz(); ++j) {
        C v = dest(0, j).value();
        if (k == 0)
          out.x[j].v = v;
        else if (k == 1)
          out.x[j].d = v;
        else
          out.x[j].dd = v;
      }
    }
    return out;
  }
  ReducedField hn(ReducedField x, int n) const {
    for (int k = 0; k < std::abs(n); ++k)
      x = h(x, n > 0);
    return x;
  }
  ReducedField dr(ReducedField x) const {
    x.p++;
    x.q++;
    for (auto &v : x.x)
      v = {v.d, v.dd, 0};
    return x;
  }
  ReducedField rho(R r) const {
    ReducedField x = field(1, 1, 0, 0);
    for (size_t j = 0; j < diff.Nz(); ++j) {
      C v = -C(1) / (r - imag_unit * a * diff.lgl_nodes()[j]);
      x.x[j] = {v, v * v, C(2) * v * v * v};
    }
    return x;
  }
  ReducedField invrho(R r) const {
    ReducedField x = field(-1, -1, 0, 0);
    for (size_t j = 0; j < diff.Nz(); ++j)
      x.x[j] = {-r + imag_unit * a * diff.lgl_nodes()[j], -1, 0};
    return x;
  }
  ReducedField tauH() const {
    return field(-1, -3, 0, 0, -imag_unit * a / std::sqrt(R(2)));
  }
  ReducedField omegaH() const {
    ReducedField x = field(-1, -1, 0, 0);
    for (size_t j = 0; j < diff.Nz(); ++j)
      x.x[j].v = -C(2) * imag_unit * a * diff.lgl_nodes()[j];
    return x;
  }
  ReducedField tau(R r) const { return mul(mul(rho(r), conj(rho(r))), tauH()); }
  ReducedField tauprimebar(R r) const {
    return C(-1) * mul(mul(conj(rho(r)), conj(rho(r))), tauH());
  }
  ReducedField eth(ReducedField x, R r) const {
    return weights(mul(conj(rho(r)), h(x)) -
                       C(x.q) *
                           mul(mul(mul(conj(rho(r)), conj(rho(r))), tauH()), x),
                   x.p + 1, x.q - 1);
  }
  // Terminal operation: only its value is consumed by the Lie components.
  ReducedField thornprime(ReducedField x, R r) const {
    ReducedField out = field(x.p - 1, x.q - 1, x.m, x.w);
    R delta = r * r - 2 * M * r + a * a;
    for (size_t j = 0; j < diff.Nz(); ++j) {
      R z = diff.lgl_nodes()[j], sig = r * r + a * a * z * z;
      C ro = -C(1) / (r - imag_unit * a * z),
        mu = std::conj(ro) * ro * ro * delta / R(2),
        gamma = mu + ro * std::conj(ro) * (r - M) / R(2);
      out.x[j].v = (-imag_unit * x.w * (r * r + a * a) / sig +
                    imag_unit * R(x.m) * a / sig - R(x.p) * gamma -
                    R(x.q) * std::conj(gamma)) *
                       x.x[j].v -
                   delta / (2 * sig) * x.x[j].d;
    }
    return out;
  }
  ReducedField rhoprime(R r) const {
    ReducedField x = field(-1, -1, 0, 0);
    for (size_t j = 0; j < diff.Nz(); ++j) {
      C ro = rho(r).x[j].v;
      x.x[j].v = -std::conj(ro) * ro * ro * (r * r - 2 * M * r + a * a) / R(2);
    }
    return x;
  }
  // One half of the paper's Hertz potential: Sdag(Phi)+conjugate,
  // matching BondiHeldMetricReconstruction's established normalization.
  ReducedField potential(ReducedField fb, R r) const {
    C P = -imag_unit * fb.w;
    ReducedField thb = conj(tauH());
    ReducedField p3 = -P * fb;
    ReducedField p2 = C(-1) * h(h(fb), false);
    ReducedField p1 = (-C(1) / (C(2) * P)) * hn(hn(fb, 2), -2) -
                      C(0.5) * (C(4) * mul(thb, h(fb)) + C(3 * M) * fb);
    ReducedField p0 =
        (-C(1) / (C(6) * P * P)) * hn(hn(fb, 3), -3) -
        (C(1) / (C(2) * P)) * (C(2) * mul(thb, h(hn(fb, 2), false)) +
                               C(M) * h(h(fb), false) - C(M) * fb);
    ReducedField ri = invrho(r);
    return weights(p0 + mul(ri, p1) + mul(mul(ri, ri), p2) +
                       mul(mul(mul(ri, ri), ri), p3),
                   -4, 0);
  }
  std::array<ReducedField, 3> primitive(ReducedField phi, R r) const {
    ReducedField ro = rho(r), rb = conj(ro), t = tau(r);
    ReducedField u = eth(phi, r) + C(3) * mul(t, phi),
                 v = dr(phi) + C(3) * mul(ro, phi);
    return {C(-1) * (eth(u, r) - mul(t, u)),
            C(-0.5) * (dr(u) - mul(ro, u) + mul(rb, u) + eth(v, r) - mul(t, v) +
                       mul(tauprimebar(r), v)),
            C(-1) * (dr(v) - mul(ro, v))};
  }
  std::array<ReducedField, 3> lie(ReducedField f, ReducedField fb, R r) const {
    C P = -imag_unit * f.w;
    ReducedField zm = h(f, false), zmb = h(fb);
    ReducedField zl = (-C(1) / (C(2) * P)) * (h(zm, false) + h(zmb));
    // Use the explicit residual-gauge construction, Hollands--Toomani
    // (J.1d). The tau terms matter in Kerr: omitting them leaves a
    // constant in h_nn.
    ReducedField zn = C(0.5) * (h(h(zl), false) + h(h(zl, false)) + zl) -
                      C(0.5) * mul(omegaH(), h(zm, false) - h(zmb)) +
                      mul(conj(tauH()), zm) + mul(tauH(), zmb);
    ReducedField ro = rho(r), rb = conj(ro), th = tauH(), thb = conj(th),
                 om = omegaH();
    ReducedField xm =
        weights(mul(conj(invrho(r)), zm) + mul(mul(ro, th), zl) - h(zl), 1, -1);
    ReducedField xn = weights(
        zn + C(0.5) * P * mul(invrho(r) + conj(invrho(r)), zl) +
            mul(mul(mul(mul(ro, rb), th), thb), zl) +
            C(0.5 * M) * mul(ro + rb, zl) + mul(mul(mul(ro, thb), om), zm) -
            mul(mul(mul(rb, th), om), zmb) - mul(mul(ro, thb), h(zl)) -
            mul(mul(rb, th), h(zl, false)),
        -1, -1);
    return {C(2) * thornprime(xn, r),
            eth(xn, r) + mul(tau(r), xn) + thornprime(xm, r) +
                mul(rhoprime(r), xm),
            C(2) * eth(xm, r)};
  }
  std::array<ReducedField, 3> rec(ReducedField f, ReducedField fb, R r) const {
    auto pos = primitive(potential(fb, r), r),
         neg = primitive(potential(conj(f), r), r);
    return {pos[0] + conj(neg[0]), conj(neg[1]), conj(neg[2])};
  }
  R norm(const ReducedField &x) const {
    R ans = 0;
    for (size_t j = 0; j < diff.Nz(); ++j)
      ans =
          std::max(ans, R(std::abs(x.x[j].v) *
                          Harm::pole_factor(x.m, x.s(), diff.lgl_nodes()[j])));
    return ans;
  }
};

using ModeCoefficients = std::map<int, C>;

ModeCoefficients held_eth_matrix_action(const ModeCoefficients &f, int spin,
                                         int m, R aw, bool up,
                                         int ell_min, int ell_max) {
  ModeCoefficients out;
  const auto coeff = [&](int ell) {
    auto it = f.find(ell);
    return it == f.end() ? C(0) : it->second;
  };
  for (int ell = ell_min; ell <= ell_max; ++ell) {
    const R l = R(ell);
    if (up) {
      const R diag = (R(1) - aw * R(m) / (l * (l + R(1)))) *
                     std::sqrt((l - R(spin)) * (l + R(spin) + R(1)) / R(2));
      const R plus = aw / (l + R(1)) *
                     std::sqrt((l - R(m) + R(1)) * (l + R(m) + R(1)) *
                               (l - R(spin)) * (l - R(spin) + R(1)) /
                               (R(2) * (R(2) * l + R(1)) * (R(2) * l + R(3))));
      const R minus = aw / l *
                      std::sqrt((l - R(m)) * (l + R(m)) *
                                (l + R(spin)) * (l + R(spin) + R(1)) /
                                (R(2) * (R(2) * l - R(1)) * (R(2) * l + R(1))));
      out[ell] = diag * coeff(ell) + plus * coeff(ell + 1) - minus * coeff(ell - 1);
    } else {
      const R diag = (R(-1) + aw * R(m) / (l * (l + R(1)))) *
                     std::sqrt((l + R(spin)) * (l - R(spin) + R(1)) / R(2));
      const R plus = aw / (l + R(1)) *
                     std::sqrt((l - R(m) + R(1)) * (l + R(m) + R(1)) *
                               (l + R(spin)) * (l + R(spin) + R(1)) /
                               (R(2) * (R(2) * l + R(1)) * (R(2) * l + R(3))));
      const R minus = aw / l *
                      std::sqrt((l - R(m)) * (l + R(m)) *
                                (l - R(spin)) * (l - R(spin) + R(1)) /
                                (R(2) * (R(2) * l - R(1)) * (R(2) * l + R(1))));
      out[ell] = diag * coeff(ell) + plus * coeff(ell + 1) - minus * coeff(ell - 1);
    }
  }
  return out;
}

ModeCoefficients omega_cos_matrix_action(const ModeCoefficients &f, int spin,
                                         int m, R a, int ell_min,
                                         int ell_max) {
  // The production Held scalar is Omega_H = -2 i a cos(theta).  Expressing
  // it through the same C_cos matrix used by the Phi Horner path gives a
  // convention-independent check of the field-level multiplication.
  ModeCoefficients out;
  const auto coeff = [&](int ell) {
    auto it = f.find(ell);
    return it == f.end() ? C(0) : it->second;
  };
  for (int ell = ell_min; ell <= ell_max; ++ell) {
    const R l = R(ell);
    const R diag = R(spin * m) / (l * (l + R(1)));
    const R plus = ghz::asymptotic::cos_theta_lower_coefficient(ell + 1, spin, m);
    const R minus = ghz::asymptotic::cos_theta_lower_coefficient(ell, spin, m);
    out[ell] = C(0, 2) * a * diag * coeff(ell) - C(0, 2) * a * plus * coeff(ell + 1) -
               C(0, 2) * a * minus * coeff(ell - 1);
  }
  return out;
}

ModeCoefficients tau_matrix_action(const ModeCoefficients &f, int spin, int m,
                                   R a, int ell_min, int ell_max) {
  // Matrix form from the adjusted operator block.  This is the coefficient
  // action of the reduced Held scalar tauH = -i a/sqrt(2), with output spin
  // s+1.  Every term, including the off-diagonal terms, carries a.
  ModeCoefficients out;
  const auto coeff = [&](int ell) {
    auto it = f.find(ell);
    return it == f.end() ? C(0) : it->second;
  };
  for (int ell = ell_min; ell <= ell_max; ++ell) {
    const R l = R(ell);
    const R diag = R(m) / (l * (l + R(1))) *
                   std::sqrt((l - R(spin)) * (l + R(spin) + R(1)) / R(2));
    const R plus = R(1) / (l + R(1)) *
                   std::sqrt((l - R(m) + R(1)) * (l + R(m) + R(1)) *
                             (l - R(spin)) * (l - R(spin) + R(1)) /
                             (R(2) * (R(2) * l + R(1)) * (R(2) * l + R(3))));
    const R minus = R(1) / l *
                    std::sqrt((l - R(m)) * (l + R(m)) *
                              (l + R(spin)) * (l + R(spin) + R(1)) /
                              (R(2) * (R(2) * l - R(1)) * (R(2) * l + R(1))));
    out[ell] = C(0, -1) * a * diag * coeff(ell) +
               C(0, 1) * a * plus * coeff(ell + 1) -
               C(0, 1) * a * minus * coeff(ell - 1);
  }
  return out;
}

ReducedField make_mode_field(const KerrDiagnostic &c, int p, int q, int m,
                             R w, const ModeCoefficients &coeffs,
                             int ell_min, int ell_max) {
  ReducedField out = c.field(p, q, m, w);
  const int spin = (p - q) / 2;
  for (size_t j = 0; j < out.x.size(); ++j) {
    for (int ell = ell_min; ell <= ell_max; ++ell) {
      auto it = coeffs.find(ell);
      if (it != coeffs.end())
        out.x[j].v += it->second * Harm::reduced_spin_weighted_spherical_harmonic(
            spin, ell, m, c.diff.lgl_nodes()[j]);
    }
  }
  return out;
}

R max_reduced_error(const ReducedField &x, const ReducedField &y) {
  if (x.x.size() != y.x.size() || x.p != y.p || x.q != y.q || x.m != y.m)
    throw std::runtime_error("operator consistency field mismatch");
  R err = 0;
  for (size_t j = 0; j < x.x.size(); ++j)
    err = std::max(err, std::abs(x.x[j].v - y.x[j].v));
  return err;
}

void check_held_operator_matrix(const KerrDiagnostic &c, int m, R w) {
  const int ell_min = 2, ell_max = 8, spin = -2;
  ModeCoefficients input{{2, C(0.31, -0.17)}, {4, C(-0.21, 0.44)},
                          {7, C(0.13, 0.09)}};
  const ReducedField x = make_mode_field(c, 0, 4, m, w, input, ell_min, ell_max);
  const R aw = c.a * w;
  for (bool up : {true, false}) {
    const ReducedField direct = c.h(x, up);
    const ModeCoefficients output =
        held_eth_matrix_action(input, spin, m, aw, up, ell_min, ell_max);
    const int output_spin = up ? spin + 1 : spin - 1;
    const int output_ell_min = std::max(ell_min, std::abs(output_spin));
    const ReducedField matrix = make_mode_field(
        c, up ? 0 : -2, up ? 2 : 4, m, w, output, output_ell_min, ell_max);
    const R error = max_reduced_error(direct, matrix);
    std::cout << "Held matrix consistency eth" << (up ? "" : "bar")
              << " error=" << std::setprecision(8) << error << "\n";
    if (!std::isfinite(error) || error > R(5e-8))
      throw std::runtime_error("Held matrix operator consistency failed");
  }
  const ReducedField omega_direct = c.mul(c.omegaH(), x);
  const ReducedField omega_matrix = make_mode_field(
      c, -1, 3, m, w, omega_cos_matrix_action(input, spin, m, c.a, ell_min,
                                               ell_max),
      ell_min, ell_max);
  const R omega_error = max_reduced_error(omega_direct, omega_matrix);
  std::cout << "Held matrix consistency OmegaH error=" << std::setprecision(8)
            << omega_error << "\n";
  if (!std::isfinite(omega_error) || omega_error > R(5e-8))
    throw std::runtime_error("Held Omega scalar consistency failed");
  const ReducedField tau_direct = c.mul(c.tauH(), x);
  const ReducedField tau_matrix = make_mode_field(
      c, -1, 1, m, w, tau_matrix_action(input, spin, m, c.a, ell_min, ell_max),
      ell_min, ell_max);
  const R tau_error = max_reduced_error(tau_direct, tau_matrix);
  std::cout << "Held matrix consistency tauH error=" << std::setprecision(8)
            << tau_error << "\n";
  if (!std::isfinite(tau_error) || tau_error > R(5e-8))
    throw std::runtime_error("Held tau scalar consistency failed");
}

struct Row {
  R a, w;
  int m, L;
  C q;
};
std::vector<Row> read(const char *path) {
  std::ifstream in(path);
  if (!in)
    throw std::runtime_error("missing input");
  std::vector<Row> rows;
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#' || line[0] == 'a')
      continue;
    std::replace(line.begin(), line.end(), ',', ' ');
    std::istringstream s(line);
    Row x;
    R r, re, im;
    if (!(s >> x.a >> r >> x.m >> x.w >> x.L >> re >> im))
      throw std::runtime_error("bad row");
    x.q = {re, im};
    rows.push_back(x);
  }
  if (rows.empty())
    throw std::runtime_error("empty input");
  return rows;
}
int main(int argc, char **argv) {
  try {
    if (argc < 3)
      throw std::runtime_error(
          "usage: bondi_kerr_m_mode_check input.csv output.csv [Nz=17]");
    auto rows = read(argv[1]);
    R a = rows.front().a, w = rows.front().w;
    int m = 2;
    size_t nz = argc > 3 ? std::stoul(argv[3]) : 17;
    if (nz < 9 || nz % 2 == 0)
      throw std::runtime_error("Nz must be odd and >=9 (equator sample)");
    spectral::SpectralDiffer diff(nz, 4);
    KerrParams kp(1, a);
    KerrMetric metric(kp);
    CoordinateHelper coords(metric);
    KinnersleyTetrad<OutgoingCoords> tetrad(metric, coords);
    OutgoingCoords point(0, 10, 0, 0);
    auto held =
        ghp::build_held_fields_vectorized(tetrad, diff.lgl_nodes(), point)
            .fields;
    ghz::numeric::AffineMap1D map(ghz::numeric::Domain{2, 20});
    KinnersleyHeldOperators<OutgoingCoords, spectral::SpectralDiffer> op(
        diff, held, tetrad, map);
    KerrDiagnostic c(diff, op, a);
    check_held_operator_matrix(c, 2, w);
    ReducedField q = c.field(0, 4, m, w), qm = c.field(0, 4, -m, -w);
    for (auto row : rows) {
      ReducedField &out = row.m == m ? q : qm;
      for (size_t j = 0; j < nz; ++j)
        out.x[j].v += row.q * Harm::reduced_spin_weighted_spherical_harmonic(
                                  -2, row.L, row.m, diff.lgl_nodes()[j]);
    }
    // Direct coupled fourth-order system: avoids squaring the angular condition
    // number through BA, and enforces both paired psi4 relations
    // simultaneously.
    ReducedField qb = conj(qm), f = c.field(4, 0, m, w),
                 fb = c.field(0, 4, m, w);
    Eigen::MatrixXcd A = Eigen::MatrixXcd::Zero(2 * nz, 2 * nz);
    Eigen::VectorXcd rhs(2 * nz);
    C coupling = C(3) * imag_unit * w;
    for (size_t j = 0; j < nz; ++j) {
      ReducedField e = f;
      e.x[j].v = 1;
      ReducedField down = c.hn(e, -4);
      e = fb;
      e.x[j].v = 1;
      ReducedField up = c.hn(e, 4);
      for (size_t k = 0; k < nz; ++k) {
        A(k, j) = down.x[k].v;
        A(nz + k, nz + j) = up.x[k].v;
      }
      A(j, nz + j) = coupling;
      A(nz + j, j) = coupling;
      rhs(j) = C(2) * imag_unit / w * q.x[j].v;
      rhs(nz + j) = C(2) * imag_unit / w * qb.x[j].v;
    }
    Eigen::VectorXd rs(2 * nz), cs(2 * nz);
    for (size_t j = 0; j < 2 * nz; ++j)
      rs(j) = 1 / A.row(j).cwiseAbs().maxCoeff();
    Eigen::MatrixXcd E = rs.asDiagonal() * A;
    for (size_t j = 0; j < 2 * nz; ++j)
      cs(j) = 1 / E.col(j).cwiseAbs().maxCoeff();
    E = E * cs.asDiagonal();
    Eigen::VectorXcd b = rs.asDiagonal() * rhs;
    auto lu = E.fullPivLu();
    if (!lu.isInvertible())
      throw std::runtime_error("singular coupled solve");
    Eigen::VectorXcd y = lu.solve(b);
    for (int k = 0; k < 3; ++k)
      y += lu.solve(b - E * y);
    Eigen::VectorXcd sol = cs.asDiagonal() * y;
    for (size_t j = 0; j < nz; ++j) {
      f.x[j].v = sol(j);
      fb.x[j].v = sol(nz + j);
    }
    R residual =
        (A * sol - rhs).cwiseAbs().maxCoeff() / rhs.cwiseAbs().maxCoeff();
    std::cout << std::setprecision(12) << "a=" << a << " Nz=" << nz
              << " omega=" << w << " coupled residual=" << residual << "\n";
    std::cout << "f(z=0)=" << f.x[nz / 2].v << " fbar(z=0)=" << fb.x[nz / 2].v
              << "\n";
    // Schwarzschild oracle checks both pieces separately, before cancellations.
    if (a == 0) {
      V fv(nz), bv(nz);
      for (size_t j = 0; j < nz; ++j) {
        fv[j] = f.x[j].v;
        bv[j] = fb.x[j].v;
      }
      ghz::asymptotic::BondiHeldMetricReconstruction old(diff, op, 1);
      auto oracle = old.reconstruct(fv, bv, m, w);
      std::array<ghz::asymptotic::MetricGridParts, 3> parts{
          oracle.nn, oracle.nm, oracle.mm};
      R worst = 0;
      for (R r : {R(20), R(100), R(500)}) {
        auto rec = c.rec(f, fb, r), lie = c.lie(f, fb, r);
        for (int k = 0; k < 3; ++k)
          for (int t = 0; t < 2; ++t) {
            ReducedField expected = t ? lie[k] : rec[k];
            auto &ser = t ? parts[k].lie_zeta : parts[k].reconstructed;
            for (size_t j = 0; j < nz; ++j)
              expected.x[j].v =
                  ser.coefficient[0][j] * r + ser.coefficient[1][j] +
                  ser.coefficient[2][j] / r + ser.coefficient[3][j] / (r * r);
            R err = c.norm((t ? lie[k] : rec[k]) - expected) / c.norm(expected);
            if (r == 20)
              std::cout << "oracle component=" << k << " lie=" << t
                        << " error=" << err << "\n";
            worst = std::max(worst, err);
          }
      }
      std::cout << "Schwarzschild separate-piece relative error=" << worst
                << "\n";
      if (worst > 1e-6)
        throw std::runtime_error("Schwarzschild oracle failed");
    }
    std::ofstream out(argv[2]);
    if (!out)
      throw std::runtime_error("cannot open radial output");
    out << std::setprecision(17)
        << "a,Nz,m,omega,r,component,reconstructed,lie,adjusted,equator_"
           "adjusted_real,equator_adjusted_imag\n";
    const char *names[] = {"nn", "nm", "mm"};
    std::array<R, 3> first{}, last{};
    for (int j = 0; j <= 60; ++j) {
      R r = 20 * std::pow(R(100), R(j) / 40);
      auto rec = c.rec(f, fb, r), lie = c.lie(f, fb, r);
      for (int k = 0; k < 3; ++k) {
        R n = c.norm(rec[k] + lie[k]);
        out << a << ',' << nz << ',' << m << ',' << w << ',' << r << ','
            << names[k] << ',' << c.norm(rec[k]) << ',' << c.norm(lie[k]) << ','
            << n << ',' << std::real((rec[k] + lie[k]).x[nz / 2].v) << ','
            << std::imag((rec[k] + lie[k]).x[nz / 2].v) << '\n';
        if (j == 20)
          first[k] = n;
        if (j == 40)
          last[k] = n;
      }
    }
    for (int k = 0; k < 3; ++k) {
      R slope = std::log(last[k] / first[k]) / std::log(R(10));
      std::cout << names[k] << " slope(200,2000)=" << slope << "\n";
      if (!std::isfinite(slope) || std::abs(slope + 1) > 0.03)
        throw std::runtime_error("radial decay check failed");
    }
    if (!std::isfinite(residual) || residual > 1e-7)
      throw std::runtime_error("seed residual too large");
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
