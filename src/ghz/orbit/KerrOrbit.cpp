//
// Created by Peter Zimmerman on 31.10.25.
//
#include "ghz/geom/KerrMetric.hpp"
#include "ghz/orbit/KerrOrbit.hpp"
#include "ghz/orbit/Splines.hpp"
#include "ghz/geom/Coords.hpp"
#include <fftw3.h>
#include <fstream>
#include <iomanip>
#ifdef _OPENMP
#include <omp.h>
#endif

// make math::sqr available
using math::sqr;

/** @brief Initializes the KerrBoundOrbit object by computing frequencies, sampling data,
* and building splines for the orbital parameters.
 * @details The initialization process involves several steps:
 * 1. Compute the orbital frequencies (f_t, f_phi, f_r, f_z) using the compute_frequencies() method. \n
 * 2. Sample the frequencies on a grid defined by Nr_ and Nz_ using sample_frequencies_for_fft(). \n
 * 3. Compute the q grid values from the sampled frequencies with compute_q_grids_from_samples(). \n
 * 4. Sample the time and azimuthal angle functions for FFT using sample_T_and_Phi_for_fft(). \n
 * 5. Re-sample the frequencies to ensure consistency. (why?) \n
 * 6. Compute the oscillatory components (Deltas) using compute_Deltas(). \n
 * 7. Finally, build periodic splines for the orbital parameters with build_splines(). \n
 *
 * */

void orbit::KerrBoundOrbit::initialize_orbit() {
    // workflow to get the orbital data
    compute_frequencies();
    sample_frequencies_for_fft(Nr_, Nz_);
    compute_q_grids_from_samples(Nr_, Nz_, qr_vals, qz_vals);
    compute_Deltas(Nr_, Nz_);
    build_splines();
}

//  evaluate orbit at a given Mino time lambda \lambda = \int d\tau/\Sigma,
//  where \tau is proper time and \Sigma = r^2 + a^2 z^2
BLCoords orbit::KerrBoundOrbit::eval_at_Mino(const Real& lambda) const {
    // linear torus angles representing mean growth of the phases from the turning points
    Real q_r = torus_freqs_.Ups_r * lambda;
    Real q_z = torus_freqs_.Ups_z * lambda;

    // interpolate instantaneous phases which fluctuate around the linear growth
    // implementation of Eq. (250) Pound and Wardell
    Real psi_r = q_r - interp_dpsi_r.eval(0.0) + interp_dpsi_r.eval(q_r);
    Real psi_z = q_z - interp_dpsi_z.eval(0.0) + interp_dpsi_z.eval(q_z);

    Real t = torus_freqs_.Ups_t * lambda
             - interp_t_r.eval(0.0) - interp_t_z.eval(0.0)
             + interp_t_r.eval(q_r) + interp_t_z.eval(q_z);

    Real phi = torus_freqs_.Ups_phi * lambda
               - interp_phi_r.eval(0.0) - interp_phi_z.eval(0.0) // spline for the oscillatory correction
               + interp_phi_r.eval(q_r) + interp_phi_z.eval(q_z);

    // return the BL coordinates evaluated at this Mino time using the interpolated phases and mean frequencies
    return BLCoords{
            t,
            p_*M_/(Real(1.0) + e_*cos(psi_r)),
            zmax_*cos(psi_z),
            phi
    };
}
/*
 * The set_constants_of_motion() method computes the specific energy (E),
 * specific angular momentum (Lz), and Carter constant (Q) for a bound orbit in Kerr spacetime,
 * based on the orbital parameters ra, rp, zmax, etc.
 */
void orbit::KerrBoundOrbit::set_constants_of_motion() {

    auto dets = computeDeterminants_(rp_, ra_);

    // see pound and wardell Eqs.(222)-(228) https://arxiv.org/pdf/2101.04592.pdf
    Real E2numer = two*dets.dg*dets.gh - dets.dh*dets.hf -two*chi_ * sqrt(
            sqr(dets.dg*dets.gh) +
            dets.hd*dets.dg*dets.gh*dets.hf + dets.hd*dets.dh*dets.hg*dets.gf );
    Real E2denom = sqr(dets.fh) + Real(4.0)*dets.fg*dets.gh;

    // specific Energy (per unit rest mass) E = -u_t = -p_t/mu
    E2_ = E2numer/E2denom;
    assert(E2_ > 0.0L && "Error: Non-physical parameters lead to complex energy!");
    E_ = sqrt(E2_);   // Energy

    assert(sqr(g_(rp_)/h_(rp_))*E2_ + (f_(rp_)*E2_-d_(rp_))/h_(rp_) >= 0.0L &&
           "Error: Non-physical parameters lead to complex angular momentum!");
    // specific angular momentum Lz = u_phi = p_phi/mu
    Lz_ = -g_(rp_)*gKerr.M()*E_/h_(rp_) + gKerr.M()*chi_ *
            sqrt(
                    sqr(g_(rp_)/h_(rp_))*E2_ + (f_(rp_)*E2_-d_(rp_))/h_(rp_)
            );
    gamma_ = one-E2_;
    // Carter constant Q  (see 6.6.1 of Pound and Wardell https://arxiv.org/pdf/2101.04592)
    Q_ = (zmax_*zmax_) * (
            sqr(gKerr.a())*gamma_ + sqr(Lz_/cos(inc_))
    );

    // auxiliary derived constants needed for root finding of R(r)=0 and frequency computations
    alpha_ = two*gKerr.M()/gamma_ - (rp_ + ra_);
    beta_  = sqr(gKerr.a())*Q_/(gamma_*rp_*ra_);

    assert(sqr(alpha_)-Real(4.0)*beta_ > 0.0L &&
           "Error: Non-physical parameters lead to complex roots!");
}

void orbit::KerrBoundOrbit::compute_torus_frequencies() {
    //
    const Real one = teuk::one;
    const Real two = teuk::two;
    const Real four = teuk::two*teuk::two;


    Real kr = (ra_-rp_)/(ra_-r3_)*(r3_-r4_)/(rp_-r4_);
    Real kz = sqr(zmax_/z1_);
    Real M = gKerr.M();
    Real a = gKerr.a();
    Real rplus = gKerr.r_plus();
    Real rminus = gKerr.r_minus();

    std::cout << "ra = " << ra_ << ", rp = " << rp_ << ", r3 = " << r3_ << ", r4 = " << r4_ << std::endl;
    std::cout << "kr = " << kr << ", kz = " << kz << std::endl;


    mino_torus_freqs.Ups_r = Real(M_PI)*sqrt(gamma_*(ra_-r3_)*(rp_-r4_)) /
                 ( two*EllipticIntegrals::computeFirstKind(kr));
    mino_torus_freqs.Ups_z = Real(M_PI)*z1_*sqrt(sqr(gKerr.a()*gamma_)) /
                 (two*EllipticIntegrals::computeFirstKind(kz));

    mino_torus_freqs.Ups_t = E_ /two * (
            r3_*(ra_+rp_+r3_)-ra_*rp_+(ra_+rp_+r3_+r4_)* Fr_(kr)+(ra_-r3_)*(rp_-r4_)* G_(kr) )
                 + two*M/(rplus-rminus) * (
            ((four*M*M*E_-a*Lz_)*rplus - two*M*a*a*E_)/(r3_-rplus)*(
                   one - F_pm_(rplus,kr)/(rp_-rplus))
            - ((four*M*M*E_-a*Lz_)*rminus - two*M*a*a*E_)/(r3_-rminus)*(
                one- F_pm_(rminus,kr)/(rp_-rminus))
    )
                 + four*M*M*E_
                 + E_*Q_*(one-G_(kz))/(gamma_*zmax_*zmax_)
                 + two*M*E_*(r3_- Fr_(kr));
    mino_torus_freqs.Ups_phi = a/(rplus-rminus) * (
            (two*M*E_*rplus-a*Lz_)/(r3_-rplus)*(one- F_pm_(rplus,kr)/(rp_-rplus))
            -(two*M*E_*rminus-a*Lz_)/(r3_-rminus)*(one- F_pm_(rminus,kr)/(rp_-rminus))
    ) + Lz_*EllipticIntegrals::computeThirdKind(zmax_*zmax_,kz)/EllipticIntegrals::computeFirstKind(kz);

    // set the torus frequencies with the Mino time frequencies
    torus_freqs_.Ups_t = mino_torus_freqs.Ups_t;
    torus_freqs_.Ups_r = mino_torus_freqs.Ups_r;
    torus_freqs_.Ups_z = mino_torus_freqs.Ups_z;
    torus_freqs_.Ups_phi = mino_torus_freqs.Ups_phi;

    // set Boyer-Lindquist mean frequencies at parameter time t_BL
    torus_freqs_.Omega_t = one;
    torus_freqs_.Omega_r = mino_torus_freqs.Ups_r/mino_torus_freqs.Ups_t;
    torus_freqs_.Omega_z = mino_torus_freqs.Ups_z/mino_torus_freqs.Ups_t;
    torus_freqs_.Omega_phi = mino_torus_freqs.Ups_phi/mino_torus_freqs.Ups_t;

} // compute_torus_frequencies
/**
 * get_T_r() returns T_r(r) (209) of Pound and Wardell
 * https://arxiv.org/pdf/2101.04592
 * @return T_r(r)
 */
Real orbit::KerrBoundOrbit::get_T_r()  const {
    Real r = p_*gKerr.M()/(one+e_*cos(phases_.psi_r));
    Real P = E_*(r*r+a_*a_)-a_*Lz_;
    //Real R = P*P - gKerr.Delta(r)*(r*r+sqr(a_*E_-Lz_)+Q_);
    return (r*r+a_*a_)/gKerr.Delta(r) * P;
}
/**
 * get_T_z() returns T_z(z) (210) of Pound and Wardell
 * https://arxiv.org/pdf/2101.04592
 * @return T_z(z)
 */
Real orbit::KerrBoundOrbit::get_T_z()  const {
    Real z = zmax_ * cos(phases_.psi_z);  // Kepler parametrize z(\psi_z)
    return -a_*a_ * E_ * (one-z*z);
}
/**
 * get_Phi_r() returns Phi_r(r) (211) of Pound and Wardell
 * https://arxiv.org/pdf/2101.04592
 * @return Phi_r(r)
 */
Real orbit::KerrBoundOrbit::get_Phi_r() const {
    Real r = p_*gKerr.M()/(one+e_*cos(phases_.psi_r));   // Kepler parametrize r(\psi_r)
    Real P = E_*(r*r+a_*a_)-a_*Lz_;
    return a_*P/gKerr.Delta(r);

}
/**
 * get_Phi_z() returns Phi_z(z) (212) of Pound and Wardell
 * https://arxiv.org/pdf/2101.04592
 * @return Phi_z(z)
 */
Real orbit::KerrBoundOrbit::get_Phi_z()  const {
    Real z = zmax_ * cos(phases_.psi_z); // Kepler parametrize z(\psi_z)
    return Lz_/(one-z*z);
}
/**
 * compute_frequencies
 * @brief Computes and sets the instantaneous frequencies d psi_alpha/d lambda = f_alpha, alpha=0,..3
 * in Mino time $lambda = int dtau/Sigma$ using the formulas in https://arxiv.org/pdf/0906.1420
 * and https://arxiv.org/pdf/2101.04592 (Eqs. 205, 206, 216, 217)
 *
 * - sets the private Real frequency variables freqs_.f_alpha
 */
void orbit::KerrBoundOrbit::compute_frequencies() {
    Real p3 = r3_*(one-e_)/M_;
    Real p4 = r4_*(one+e_)/M_;
    freqs_.f_t = get_T_r() + get_T_z() + a_*Lz_;  // Mino f_t = dt/d\lambda = Sigma dt/d\tau
    freqs_.f_phi = get_Phi_r() + get_Phi_z() - a_*Lz_;  // Mino f_phi = dphi/d\lambda = Sigma dphi/d\tau

    Real term1 = p_ - p3 - e_*(p_+p3*cos(phases_.psi_r));
    Real term2 = p_ - p4 + e_*(p_-p4*cos(phases_.psi_r));
    freqs_.f_r = M_ * sqrt(gamma_*term1*term2)/(one-e_*e_);  // (216) Pound and Wardell

    freqs_.f_z = a_*sqrt(
            gamma_ *
            (z1_*z1_-zmax_*zmax_*sqr(cos(phases_.psi_z)))
    );
}
/**
 * update_q_angles
 * @brief Update the mean action angles q_alpha
 *
 * - updates the private internal action angles q_alpha = Upsilon_alpha * lambda + q_alpha^0
 */
void orbit::KerrBoundOrbit::update_q_angles(orbit::Real const& lambda) {
    torus_angles_.q_r = torus_angles_.q_r + torus_freqs_.Ups_r*lambda;
    torus_angles_.q_z = torus_angles_.q_z + torus_freqs_.Ups_z*lambda;
    torus_angles_.q_t = torus_angles_.q_t + torus_freqs_.Ups_t*lambda;
    torus_angles_.q_phi = torus_angles_.q_phi + torus_freqs_.Ups_phi*lambda;
}

std::vector<Complex> orbit::KerrBoundOrbit::resample_to_uniform_q(
        const std::vector<Complex>& samples_in,
        const std::vector<Real>& q_grid) const
{
    const std::size_t N = samples_in.size();
    if (q_grid.size() != N) {
        throw std::runtime_error("resample_to_uniform_q: size mismatch");
    }

    std::vector<Real> re_in(N), im_in(N);
    for (std::size_t i = 0; i < N; ++i) {
        re_in[i] = samples_in[i].real();
        im_in[i] = samples_in[i].imag();
    }

    PeriodicSpline spline_re(q_grid, re_in, teuk::twoPi);
    PeriodicSpline spline_im(q_grid, im_in, teuk::twoPi);

    std::vector<Complex> out(N);
    for (std::size_t j = 0; j < N; ++j) {
        Real q = teuk::twoPi * Real(j) / Real(N);
        out[j] = Complex(spline_re.eval(q), spline_im.eval(q));
    }
    return out;
}

void orbit::KerrBoundOrbit::sample_frequencies_for_fft(size_t Nr, size_t Nz) {

    const Real psi_r_saved = phases_.psi_r;
    const Real psi_z_saved = phases_.psi_z;

    f_samples_.resize(Nr*Nz);

    for (size_t i=0; i<Nr; ++i) {
        Real psi_r = teuk::twoPi*i/Nr;

        for (size_t j=0; j<Nz; ++j) {
            Real psi_z = teuk::twoPi*j/Nz;

            // Set phases
            phases_.psi_r = psi_r;
            phases_.psi_z = psi_z;

            // Compute instantaneous frequencies at these phases
            compute_frequencies();

            // Flat index into 2D grid
            size_t idx = i*Nz + j; // i is the r-index and j is the z-index

            // Store full frequency set on the uniform grid
            f_samples_[idx] = {freqs_.f_t, freqs_.f_phi, freqs_.f_r, freqs_.f_z};
        }
    }
    // restore original orbit state!
    phases_.psi_r = psi_r_saved;
    phases_.psi_z = psi_z_saved;
}

// Input: f_r[i] and f_z[j] extracted from f_samples_
// Output: q_r_uniform[i] and q_z_uniform[j] ∈ [0,2π)

void orbit::KerrBoundOrbit::compute_q_from_psi(const std::vector<Real>& f, std::vector<Real>& q) {
    size_t N = f.size();

    std::vector<Real> invf(N), cum(N+1,0.0L);
    for (size_t i=0;i<N;i++)
        invf[i] = 1.0L / f[i]; //  invf tells how much affine parameter λ changes per unit ψ

    // cumulative integral cum = the cumulative λ along the orbit
    for (size_t i=0;i<N;i++){
        Real dpsi = twoPi / N; // because ψ-grid is uniform
        cum[i+1] = cum[i] + 0.5L*(invf[i] + invf[(i+1)%N])*dpsi; // crude trapezoidal rule
        // Notice the modulo (i+1)%N makes it periodic
    }
    // Normalize to get q = Upsilon * cum, rescaled to [0,2π)
    Real scale = twoPi / cum[N]; // = Upsilon / <f> automatically since cum[N] is the total "affine length" over one cycle of ψ
    q.resize(N);
    for (size_t i=0;i<N;i++)
        q[i] = scale * cum[i]; // set q[i] s.t. it increases uniformly in torus phase advance
}

void orbit::KerrBoundOrbit::compute_q_grids_from_samples(
        const size_t& Nr, const size_t& Nz,
        std::vector<Real>& q_r_uniform,
        std::vector<Real>& q_z_uniform)
{
    std::vector<Real> f_r_slice(Nr), f_z_slice(Nz);

    for (size_t i = 0; i < Nr; ++i) {
        size_t idx = i * Nz;
        f_r_slice[i] = f_samples_[idx].f_r;
    }
    for (size_t j = 0; j < Nz; ++j) {
        size_t idx = j;
        f_z_slice[j] = f_samples_[idx].f_z;
    }

    compute_q_from_psi(f_r_slice, q_of_psi_r_vals);
    compute_q_from_psi(f_z_slice, q_of_psi_z_vals);

    q_r_uniform.resize(Nr);
    q_z_uniform.resize(Nz);
    for (size_t i = 0; i < Nr; ++i) q_r_uniform[i] = teuk::twoPi * Real(i) / Real(Nr);
    for (size_t j = 0; j < Nz; ++j) q_z_uniform[j] = teuk::twoPi * Real(j) / Real(Nz);
}

// map FFT index -> integer mode k in range -(N/2)..(N/2-1)
inline int fft_index_to_mode(int idx, int N) {
    if (idx <= N/2) return idx;
    return idx - N;
}

// Assumes samples_q are real-valued data stored in Complex form,
// so negative Fourier modes are reconstructed by conjugation.
void orbit::KerrBoundOrbit::compute_Delta_and_freq_modes_SIMD(
        const std::vector<Complex>& samples_in,
        const std::vector<Real>& q_of_psi,
        const Real& Ups,
        std::vector<Real>& Delta_out,
        std::vector<Complex>& modes_out)
{
    using LD = long double;
    using LC = std::complex<long double>;

    const size_t N = samples_in.size();
    const LC Ilc(0, 1);
    const LD UpsL = LD(Ups);
    /*
        Algorithm
        ---------------------------
        for fixed j (fixed q_z)
        take samples_in[i = 0..Nr-1]    // slice along r
           use q_r_vals[i] in reconstruction
           1D FFT → modes f_{k_r}(q_z)
           end j
       then for each k_r:
           take values vs j
           use q_z_vals[j]
       1D FFT → f_{k_r,k_z}
       divide by -i (k_r Ω_r) or  - i (k_z Ω_z)
       inverse FFT back
   */
    // FFTW uses long double precision
    // unfortunately there is no multi-precision complex type in FFTW

    // 1. Resample X(psi_i) -> X(q_j) on uniform q-grid
    std::vector<Complex> samples_q = resample_to_uniform_q(samples_in, q_of_psi);

    // 2. Compute mean on uniform q-grid
    Complex mean = std::accumulate(samples_q.begin(), samples_q.end(), Complex(0.,0.)) / Complex(N);

    // 3. Choose persistent FFT plan
    fftwl_plan plan = nullptr;
    fftwl_complex* in = nullptr;
    fftwl_complex* out = nullptr;

    if (N == Nr_) {
        plan = fft_plan_r_;
        in   = fft_in_r_;
        out  = fft_out_r_;
    } else if (N == Nz_) {
        plan = fft_plan_z_;
        in   = fft_in_z_;
        out  = fft_out_z_;
    } else {
        throw std::runtime_error("compute_Delta_and_freq_modes_SIMD: size mismatch");
    }

    // 4. Fill FFT input from uniform q samples
#pragma omp parallel for default(none) shared(in,samples_q,mean,N)
    for (size_t i = 0; i < N; ++i) {
        in[i][0] = LD(samples_q[i].real() - mean.real());
        in[i][1] = LD(samples_q[i].imag() - mean.imag());
    }

    fftwl_execute(plan);

    modes_out.assign(N, Complex(0.0, 0.0));
#pragma omp parallel for default(none) shared(out,modes_out,N)
    for (size_t k = 0; k < N; ++k) {
        modes_out[k] = LC(out[k][0], out[k][1]) / LD(N);
    }

    // 5. Reconstruct Delta on the same uniform q-grid
    Delta_out.assign(N, Real(0));
    const int kmax = static_cast<int>((N + 1) / 2);
    // Intentionally omit Nyquist mode for spectral integration.

#pragma omp parallel for default(none) shared(Delta_out,modes_out,N,kmax,Ilc,UpsL)
    for (size_t i = 0; i < N; ++i) {
        const LD q = LD(Real(2*M_PI*i) / Real(N));
        LD sum_re = 0.0L;

#pragma omp simd reduction(+:sum_re)
        for (int k = 1; k < kmax; ++k) {
            LD kL = LD(k);
            LC mk(LD(modes_out[k].real()), LD(modes_out[k].imag()));

            LD ang = -kL * q;
            LD c = cosl(ang);
            LD s = sinl(ang);

            LC e_pos(c, s);
            LC e_neg(c, -s);

            LC denom(0.0L, -kL * UpsL);
            LC term = (mk * e_neg - std::conj(mk) * e_pos) / denom;
            sum_re += term.real();
        }

        Delta_out[i] = Real(sum_re);
    }
}


/**
*   \function compute_Deltas
 *  computes the oscillating pieces (zero average) of the orbital phases
 *  as inverse FFTs
 */
void orbit::KerrBoundOrbit::compute_Deltas(size_t Nr, size_t Nz)
{
    const Real two_pi = teuk::twoPi;
    // allocate slice containers with correct sizes
    std::vector<Complex> slice_r(Nr), slice_t_r(Nr), slice_phi_r(Nr);
    std::vector<Complex> slice_z(Nz), slice_t_z(Nz), slice_phi_z(Nz);

    // radial slice (j=0)
#pragma omp parallel for default(none) shared(slice_r,slice_t_r,slice_phi_r,Nr,Nz)
    for (size_t i = 0; i < Nr; ++i) {
        size_t idx = i * Nz + 0; // j=0 slice
        slice_r[i]     = Complex(f_samples_[idx].f_r, 0.0L);
        slice_t_r[i]   = Complex(f_samples_[idx].f_t, 0.0L);
        slice_phi_r[i] = Complex(f_samples_[idx].f_phi, 0.0L);
    }

    // polar slice (choose i=0)
#pragma omp parallel for default(none) shared(slice_z,slice_t_z,slice_phi_z,Nz)
    for (size_t j = 0; j < Nz; ++j) {
        size_t idx = 0 * Nz + j; // /i=0 slice
        slice_z[j]     = Complex(f_samples_[idx].f_z, 0.0L);
        slice_t_z[j]   = Complex(f_samples_[idx].f_t, 0.0L);
        slice_phi_z[j] = Complex(f_samples_[idx].f_phi, 0.0L);
    }

    // Prepare outputs
    Delta_psi_r_.resize(Nr); Delta_psi_z_.resize(Nz);
    Delta_t_r_.resize(Nr);   Delta_t_z_.resize(Nz);
    Delta_phi_r_.resize(Nr); Delta_phi_z_.resize(Nz);

    compute_Delta_and_freq_modes_SIMD(slice_r, q_of_psi_r_vals, torus_freqs_.Ups_r,
                                      Delta_psi_r_, f_modes_.f_r_modes);
    compute_Delta_and_freq_modes_SIMD(slice_t_r, q_of_psi_r_vals, torus_freqs_.Ups_r,
                                      Delta_t_r_,   f_modes_.T_r_modes);
    compute_Delta_and_freq_modes_SIMD(slice_phi_r, q_of_psi_r_vals, torus_freqs_.Ups_r,
                                      Delta_phi_r_, f_modes_.Phi_r_modes);

    compute_Delta_and_freq_modes_SIMD(slice_z, q_of_psi_z_vals, torus_freqs_.Ups_z,
                                      Delta_psi_z_, f_modes_.f_z_modes);
    compute_Delta_and_freq_modes_SIMD(slice_t_z, q_of_psi_z_vals, torus_freqs_.Ups_z,
                                      Delta_t_z_,   f_modes_.T_z_modes);
    compute_Delta_and_freq_modes_SIMD(slice_phi_z, q_of_psi_z_vals, torus_freqs_.Ups_z,
                                      Delta_phi_z_, f_modes_.Phi_z_modes);
}

void orbit::KerrBoundOrbit::export_trajectory_stream(const std::string& filename,
                                                     const Real& lambda_max,      // maximum Mino time
                                                     const Real& dlambda,         // step size in lambda
                                                     size_t output_stride)        // output every N-th step to reduce file size

{
    std::ofstream out(filename);
    if (!out) {
        throw std::runtime_error("Cannot open file for writing trajectory.");
    }

    // Header
    out << "lambda,t,r,z,phi\n";
    out << std::fixed << std::setprecision(15);

    size_t step = 0;
    for (Real lambda = 0; lambda <= lambda_max; lambda += dlambda, ++step) {
        if (step % output_stride != 0) continue;  // skip to reduce file size

        BLCoords coords = eval_at_Mino(lambda);
        Real z = std::clamp(coords.x2, Real(-1.0), Real(1.0));

        out << lambda << ","
            << coords.x0 << ","
            << coords.x1 << ","
            << z << ","
            << coords.x3 << "\n";
    }

    out.close();
}


void orbit::KerrBoundOrbit::prepare_fft_plans() {

    if (fft_plan_r_ != nullptr) {
        fftwl_destroy_plan(fft_plan_r_);
        fftwl_free(fft_in_r_);
        fftwl_free(fft_out_r_);
        fft_plan_r_ = nullptr;
        fft_in_r_ = nullptr;
        fft_out_r_ = nullptr;
    }
    if (fft_plan_z_ != nullptr) {
        fftwl_destroy_plan(fft_plan_z_);
        fftwl_free(fft_in_z_);
        fftwl_free(fft_out_z_);
        fft_plan_z_ = nullptr;
        fft_in_z_ = nullptr;
        fft_out_z_ = nullptr;
    }

    if (fft_plan_r_ == nullptr) {
        fft_in_r_  = (fftwl_complex*) fftwl_malloc(sizeof(fftwl_complex) * Nr_);
        fft_out_r_ = (fftwl_complex*) fftwl_malloc(sizeof(fftwl_complex) * Nr_);
        fft_plan_r_ = fftwl_plan_dft_1d(static_cast<int>(Nr_), fft_in_r_, fft_out_r_, FFTW_FORWARD, FFTW_ESTIMATE);
    }

    if (fft_plan_z_ == nullptr) {
        fft_in_z_  = (fftwl_complex*) fftwl_malloc(sizeof(fftwl_complex) * Nz_);
        fft_out_z_ = (fftwl_complex*) fftwl_malloc(sizeof(fftwl_complex) * Nz_);
        fft_plan_z_ = fftwl_plan_dft_1d(static_cast<int>(Nz_), fft_in_z_, fft_out_z_, FFTW_FORWARD, FFTW_ESTIMATE);
    }

}

void orbit::KerrBoundOrbit::free_fft() {
    if (fft_plan_r_) fftwl_destroy_plan(fft_plan_r_);
    if (fft_in_r_)   fftwl_free(fft_in_r_);
    if (fft_out_r_)  fftwl_free(fft_out_r_);
    if (fft_plan_z_) fftwl_destroy_plan(fft_plan_z_);
    if (fft_in_z_)   fftwl_free(fft_in_z_);
    if (fft_out_z_)  fftwl_free(fft_out_z_);
    fft_plan_r_ = nullptr; fft_in_r_ = nullptr;
    fft_out_r_ = nullptr; fft_plan_z_ = nullptr;
    fft_in_z_ = nullptr; fft_out_z_ = nullptr;
}