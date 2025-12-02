//
// Created by Peter Zimmerman on 31.10.25.
//
#include "../include/KerrMetric.hpp"
#include "../include/KerrOrbit.hpp"
#include "../include/Splines.hpp"
#include "../include/Coords.hpp"
#include <fftw3.h>
#include <fstream>
#include <iomanip>
#ifdef _OPENMP
#include <omp.h>
#endif

// make math::sqr available
using math::sqr;

void orbit::KerrBoundOrbit::prepare_fft_plans() {

    if (fft_plan_r_ != nullptr) {
        fftwl_destroy_plan(fft_plan_r_);
        fftwl_free(fft_in_r_);
        fftwl_free(fft_out_r_);
    }
    if (fft_plan_z_ != nullptr) {
        fftwl_destroy_plan(fft_plan_z_);
        fftwl_free(fft_in_z_);
        fftwl_free(fft_out_z_);
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
void orbit::KerrBoundOrbit::init() {
    // workflow to get the orbital data
    compute_torus_frequencies(); // compute Upsilons and Omega (avg frequencies in mino and BL time)
    compute_frequencies(); // compute initial frequencies f_t, f_phi, f_r, f_z where f = d\psi/d\lambda
    sample_frequencies_for_fft(Nr_, Nz_);
    compute_q_grids_from_samples(Nr_, Nz_, qr_vals, qz_vals); // fill the q grid values
    sample_T_and_Phi_for_fft(Nr_, Nz_);
    sample_frequencies_for_fft(Nr_, Nz_);
    compute_Deltas(Nr_, Nz_);

    // build periodic splines
    build_splines();
}

//  evaluate orbit at Mino time lambda
BLCoords orbit::KerrBoundOrbit::eval_at_Mino(const Real& lambda) const {
    // linear torus angles
    Real q_r = torus_freqs_.Ups_r * lambda;
    Real q_z = torus_freqs_.Ups_z * lambda;
    // total phases
    Real psi_r = interp_psi_r.eval(q_r) + interp_dpsi_r.eval(q_r);
    Real psi_z = interp_psi_z.eval(q_z) + interp_dpsi_z.eval(q_z);

    return BLCoords{torus_freqs_.Omega_t*lambda+interp_t_r.eval(q_r)+interp_t_z.eval(q_z),
                    p_*M_/(1.0+e_*cos(psi_r)),
                    acos(zmax_*cos(psi_z)),
                    torus_freqs_.Omega_phi*lambda+interp_phi_r.eval(q_r)+interp_phi_z.eval(q_z)
    };
}

void orbit::KerrBoundOrbit::set_constants_of_motion() {

    auto dets = computeDeterminants_(rp_, ra_);

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
    Lz_ = g_(rp_)*gKerr.M()*E_/h_(rp_) + gKerr.M()*chi_ *
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


    // set Boyer-Lindquist mean frequencies at parameter time t_BL
    torus_freqs_.Omega_t = one;
    torus_freqs_.Omega_r = mino_torus_freqs.Ups_r/mino_torus_freqs.Ups_t;
    torus_freqs_.Omega_z = mino_torus_freqs.Ups_z/mino_torus_freqs.Ups_t;
    torus_freqs_.Omega_phi = mino_torus_freqs.Ups_phi/mino_torus_freqs.Ups_t;

}
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
void orbit::KerrBoundOrbit::update_q_angles(orbit::Real const& mino_time_param) {
    torus_angles_.q_r = torus_angles_.q_r + torus_freqs_.Omega_r*mino_time_param;
    torus_angles_.q_z = torus_angles_.q_z + torus_freqs_.Omega_z*mino_time_param;
    torus_angles_.q_t = torus_angles_.q_t + torus_freqs_.Omega_t*mino_time_param;
    torus_angles_.q_phi = torus_angles_.q_phi + torus_freqs_.Omega_phi*mino_time_param;
}
/**
 * sample_T_and_Phi_for_fft
 * @brief generates T_r(r) and T_z(z) for f_t
 *  and Phi_r(r) and Phi_z(z) for f_phi
 * - updates the private internal vectors T_a_samples_ and Phi_a_samples_
 */
void orbit::KerrBoundOrbit::sample_T_and_Phi_for_fft(size_t Nr, size_t Nz)
{
    const Real two_pi = teuk::twoPi;
    const Real psi_r_saved = phases_.psi_r;
    const Real psi_z_saved = phases_.psi_z;

    T_r_samples_.resize(Nr);
    Phi_r_samples_.resize(Nr);
    T_z_samples_.resize(Nz);
    Phi_z_samples_.resize(Nz);

    // Sample T_r(psi_r),Phi_r(psi_r)
    for (size_t i = 0; i < Nr; ++i) {
        Real psi_r = teuk::twoPi * i / Nr;
        phases_.psi_r = psi_r;
        // leave psi_z as whatever; T_r,Phi_r depend only on psi_r
        T_r_samples_[i] = get_T_r();
        Phi_r_samples_[i] = get_Phi_r();
    }

    // Sample T_z(psi_z),Phi_z(psi_z)
    for (size_t j = 0; j < Nz; ++j) {
        Real psi_z =teuk::twoPi * j / Nz;
        phases_.psi_z = psi_z;
        // T_z depends only on psi_z
        T_z_samples_[j] = get_T_z();
        Phi_z_samples_[j] = get_Phi_z();
    }
    // restore original orbit state!
    phases_.psi_r = psi_r_saved;
    phases_.psi_z = psi_z_saved;
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

void orbit::KerrBoundOrbit::compute_q_grids_from_samples(const size_t& Nr, const size_t& Nz,
                                                         std::vector<Real>& q_r,
                                                         std::vector<Real>& q_z)
{
    const Real twoPi = teuk::twoPi;

    // --- Extract f_r(ψ_r, ψ_z=0) → f_r_slice[i] ---
    std::vector<Real> f_r_slice(Nr);
    for (size_t i = 0; i < Nr; ++i) {
        size_t idx = i * Nz + 0;   // j = 0 slice
        f_r_slice[i] = f_samples_[idx].f_r;
    }

    // --- Extract f_z(ψ_r=0, ψ_z) → f_z_slice[j] ---
    std::vector<Real> f_z_slice(Nz);
    for (size_t j = 0; j < Nz; ++j) {
        size_t idx = 0 * Nz + j;   // i = 0 slice
        f_z_slice[j] = f_samples_[idx].f_z;
    }

    // Now compute q_r(ψ_r) and q_z(ψ_z)
    compute_q_from_psi(f_r_slice, q_r);
    compute_q_from_psi(f_z_slice, q_z);
}

// map FFT index -> integer mode k in range -(N/2)..(N/2-1)
inline int fft_index_to_mode(int idx, int N) {
    if (idx <= N/2) return idx;
    return idx - N;
}

void orbit::KerrBoundOrbit::compute_Delta_and_freq_modes(
        const std::vector<Complex>& samples_in,  // f(ψ) samples
        const std::vector<Real>& q_grid,         // q[ψ] sampled at same ψ-grid points <-- pass q_r_vals OR q_z_vals here
        const Real& Omega,                              // \Upsilon_α or Ω_α (torus freq)
        std::vector<Real>& Delta_out,            // Δψ or Δt or Δφ oscillating piece (changed in place)
        std::vector<Complex>& modes_out)         // Fourier modes c_k of frequency f_α (changed in place)
{
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
    using LD = long double;
    using LC = std::complex<long double>;

    size_t N = samples_in.size();

    Delta_out.assign(N, 0);
    modes_out.assign(N, LC(0,0));
    const LC Ilc(0, 1);

    //Delta_out.assign(N, Real(0));
    //modes_out.assign(N, Complex(0,0));

    // compute mean
    Complex mean = std::accumulate(samples_in.begin(), samples_in.end(),
                                   Complex(0.,0.) ) /
                                           Complex(N);

    std::vector<Complex> osc(N); // oscillating piece of the phase values
    // mean subtraction
#pragma omp parallel for default(none) shared(osc,samples_in,mean,N)
    for(size_t i=0; i<N;i ++) osc[i] = samples_in[i] - mean;

    // Allocate FFTW buffers
    auto* in  = (fftwl_complex*)fftwl_malloc(sizeof(fftwl_complex) * N);
    auto* out = (fftwl_complex*)fftwl_malloc(sizeof(fftwl_complex) * N);

// Fill input
#pragma omp parallel for default(none) shared(in,samples_in,mean,N)
    for(size_t i=0; i<N; i++) {
        in[i][0] = LD(samples_in[i].real() - mean.real());
        in[i][1] = LD(samples_in[i].imag() - mean.imag());
    }

// Create plan
    fftwl_plan fftwlPlan = fftwl_plan_dft_1d((int)N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
// Execute FFT
    fftwl_execute(fftwlPlan);

    //std::vector<Complex> raw_modes(N); // FFT modes
    //fftwl_plan plan = fftwl_plan_dft_1d(int(N), reinterpret_cast<fftwl_complex*>(osc.data()), reinterpret_cast<fftwl_complex*>(raw_modes.data()), FFTW_FORWARD, FFTW_ESTIMATE);
    //fftwl_execute(plan);
    //fftwl_destroy_plan(plan);

    // normalize
#pragma omp parallel for default(none) shared(out,modes_out,N)
    for(size_t k=0; k<N; k++) {
        modes_out[k] = LC(out[k][0], out[k][1]) / LD(N);
    }

    //for(size_t k=0;k<N;k++) modes_out[k] = raw_modes[k];

    // construct Δ using *actual* q-grid
    int kmax = int(N/2);
#pragma omp parallel for default(none) shared(Delta_out,modes_out,q_grid,Omega,N,kmax,Ilc)
    for(size_t i=0; i<N; i++) { // over grid points
        // get q value at this grid point
        const LD& q = (LD)q_grid[i];
        LC sum(0.0L,0.0L); // start sum at zero for grid point i

        for (int k = 1; k < kmax; ++k) { // sum over modes

            LC kC = LC((LD)k, 0.0L);

            // promote FFTW mode to long double
            LC mk = LC((LD)modes_out[k].real(), (LD)modes_out[k].imag());

            sum += mk * exp(-Ilc * kC * q)
                   / (-Ilc * kC * (LD)Omega);

            sum += std::conj(mk) * exp(Ilc * kC * q)
                   / ( Ilc * kC * (LD)Omega); // add negative modes via conj symmetry
        }

        Delta_out[i] = sum.real();
    }
    // Cleanup
    fftwl_destroy_plan(fftwlPlan);
    fftwl_free(in);
    fftwl_free(out);
}

void orbit::KerrBoundOrbit::compute_Delta_and_freq_modes_SIMD(
        const std::vector<Complex>& samples_in,  // f(ψ) samples
        const std::vector<Real>& q_grid,         // q[ψ] sampled at same ψ-grid points <-- pass q_r_vals OR q_z_vals here
        const Real& Omega,                              // \Upsilon_α or Ω_α (torus freq)
        std::vector<Real>& Delta_out,            // Δψ or Δt or Δφ oscillating piece (changed in place)
        std::vector<Complex>& modes_out)         // Fourier modes c_k of frequency f_α (changed in place)
{
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
    using LD = long double;
    using LC = std::complex<long double>;

    size_t N = samples_in.size();

    Delta_out.assign(N, 0);
    modes_out.assign(N, LC(0,0));
    const LC Ilc(0, 1);
    const LD OmL = (LD)Omega;

    //Delta_out.assign(N, Real(0));
    //modes_out.assign(N, Complex(0,0));

    // compute mean
    Complex mean = std::accumulate(samples_in.begin(), samples_in.end(),
                                   Complex(0.,0.) ) /
                   Complex(N);

    // choose FFT plan based on N
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
        throw std::runtime_error("compute_Delta_and_freq_modes: input size does not match Nr or Nz!");
    }

// Fill input with oscillating piece by subtracting mean
#pragma omp parallel for default(none) shared(in,samples_in,mean,N)
    for(size_t i=0; i<N; i++) {
        in[i][0] = LD(samples_in[i].real() - mean.real());
        in[i][1] = LD(samples_in[i].imag() - mean.imag());
    }
    fftwl_execute(plan);  // execute persistent plan

    // normalize by N size
#pragma omp parallel for default(none) shared(out,modes_out,N)
    for(size_t k=0; k<N; k++) {
        modes_out[k] = LC(out[k][0], out[k][1]) / LD(N);
    }

    // construct Δ using *actual* q-grid
    int kmax = int(N/2);
#pragma omp parallel for default(none) shared(Delta_out,modes_out,q_grid,OmL,N,kmax,Ilc)
    for(size_t i=0; i<N; i++) { // over grid points
        // get q value at this grid point
        const LD& q = (LD)q_grid[i];
        LD sum_re = 0.0L; // start sum at zero for grid point i

#pragma omp simd reduction(+:sum_re)
        for (int k = 1; k < kmax; ++k) { // sum over modes

            LD kL = (LD)k;
            const LC mk(
                    (LD)modes_out[k].real(),
                    (LD)modes_out[k].imag()
            );

            const LD ang = -kL * q;   // argument of exp
            const LD cos_a = cosl(ang);
            const LD sin_a = sinl(ang);
            // e^{± i k q}
            const LC e_pos(cos_a,  sin_a);
            const LC e_neg(cos_a, -sin_a);

            const LC denom(0.0L, -kL * OmL); // -i k Ω
            const LC term = (mk * e_neg - conj(mk) * e_pos) / denom;

            sum_re += term.real();

        }

        Delta_out[i] = sum_re;
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

    // compute radial deltas and freq modes
    compute_Delta_and_freq_modes_SIMD(slice_r, qr_vals, torus_freqs_.Ups_r,
                                 Delta_psi_r_, f_modes_.f_r_modes);
    compute_Delta_and_freq_modes_SIMD(slice_t_r, qr_vals, torus_freqs_.Ups_r,
                                 Delta_t_r_,   f_modes_.T_r_modes);
    compute_Delta_and_freq_modes_SIMD(slice_phi_r,qr_vals, torus_freqs_.Omega_r,
                                 Delta_phi_r_, f_modes_.Phi_r_modes);

    // compute polar deltas and freq modes
    compute_Delta_and_freq_modes_SIMD(slice_z, qz_vals,  torus_freqs_.Omega_z,
                                 Delta_psi_z_, f_modes_.f_z_modes);
    compute_Delta_and_freq_modes_SIMD(slice_t_z, qz_vals, torus_freqs_.Omega_z,
                                 Delta_t_z_,   f_modes_.T_z_modes);
    compute_Delta_and_freq_modes_SIMD(slice_phi_z, qz_vals, torus_freqs_.Omega_z,
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
    out << "lambda,t,r,theta,phi\n";
    out << std::fixed << std::setprecision(15);

    size_t step = 0;
    for (Real lambda = 0; lambda <= lambda_max; lambda += dlambda, ++step) {
        if (step % output_stride != 0) continue;  // skip to reduce file size

        BLCoords coords = eval_at_Mino(lambda);
        Real t = torus_freqs_.Omega_t * lambda + interp_t_r.eval(torus_freqs_.Omega_r*lambda)
                 + interp_t_z.eval(torus_freqs_.Omega_z*lambda);  // optional: total t

        out << lambda << ","
            << coords.x0 << ","
            << coords.x1 << ","
            << acos(coords.x2) << ","
            << coords.x3 << "\n";
    }

    out.close();
}