//
// Created by Peter Zimmerman on 31.10.25.
//
#include "../include/KerrMetric.hpp"
#include "../include/KerrOrbit.hpp"
#include "../include/EllipticIntegrals.hpp"
#include <fftw3.h>



void ghz::KerrBoundOrbit::set_constants_of_motion() {
    auto dets = computeDeterminants_(rp_, ra_);

    const Real one = teuk::one;
    const Real two = teuk::two;

    Real E2numer = Real(2.0)*dets.dg*dets.gh - dets.dh*dets.hf
                   - Real(2.0)*chi_*math::sqrt(
            math::sqr(dets.dg*dets.gh)+dets.hd*dets.dg*dets.gh*dets.hf+dets.hd*dets.dh*dets.hg*dets.gf
    );
    Real E2denom = math::sqr(dets.fh) + Real(4.0)*dets.fg*dets.gh;
    E2_ = E2numer/E2denom;
    E_ = math::sqrt(E2_);

    Lz_ = g_(rp_)*gKerr.M()*E_/h_(rp_) + gKerr.M()*chi_*math::sqr(
            math::sqr(g_(rp_)/h_(rp_))*E2_ + (f_(rp_)*E2_-d_(rp_))/h_(rp_)
    );
    gamma_ = teuk::one-E2_;
    Q_ = math::sqrt(zmax_) * (
            math::sqr(gKerr.a())*gamma_ + math::sqr(Lz_/math::cos(inc_))
    );

    alpha_ = two*gKerr.M()/gamma_ - rp_ + ra_;
    beta_  = math::sqr(gKerr.a())*Q_/(gamma_*rp_*ra_);
}

void ghz::KerrBoundOrbit::compute_torus_frequencies() {
    //
    const Real one = teuk::one;
    const Real two = teuk::two;
    const Real four = teuk::two*teuk::two;


    Real kr = (ra_-rp_)/(ra_-r3_)*(r3_-r4_)*(rp_-r4_);
    Real kz = math::sqr(zmax_/z1_);
    Real M = gKerr.M();
    Real a = gKerr.a();
    Real rplus = gKerr.r_plus();
    Real rminus = gKerr.r_minus();


    mino_torus_freqs.Ups_r = Real(M_PI)*math::sqrt(gamma_*(ra_-r3_)*(rp_-r4_)) /
                 ( two*EllipticIntegrals::computeFirstKind(kr));
    mino_torus_freqs.Ups_z = Real(M_PI)*z1_*math::sqrt(math::sqr(gKerr.a()*gamma_)) /
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
Real ghz::KerrBoundOrbit::get_T_r()  const {
    Real r = p_*gKerr.M()/(one+e_*math::cos(phases_.psi_r));
    Real P = E_*(r*r+a_*a_)-a_*Lz_;
    //Real R = P*P - gKerr.Delta(r)*(r*r+math::sqr(a_*E_-Lz_)+Q_);
    return (r*r+a_*a_)/gKerr.Delta(r) * P;
}
/**
 * get_T_z() returns T_z(z) (210) of Pound and Wardell
 * https://arxiv.org/pdf/2101.04592
 * @return T_z(z)
 */
Real ghz::KerrBoundOrbit::get_T_z()  const {
    Real z = zmax_ * math::cos(phases_.psi_z);  // Kepler parametrize z(\psi_z)
    return -a_*a_ * E_ * (one-z*z);
}
/**
 * get_Phi_r() returns Phi_r(r) (211) of Pound and Wardell
 * https://arxiv.org/pdf/2101.04592
 * @return Phi_r(r)
 */
Real ghz::KerrBoundOrbit::get_Phi_r() const {
    Real r = p_*gKerr.M()/(one+e_*math::cos(phases_.psi_r));   // Kepler parametrize r(\psi_r)
    Real P = E_*(r*r+a_*a_)-a_*Lz_;
    return a_*P/gKerr.Delta(r);

}
/**
 * get_Phi_z() returns Phi_z(z) (212) of Pound and Wardell
 * https://arxiv.org/pdf/2101.04592
 * @return Phi_z(z)
 */
Real ghz::KerrBoundOrbit::get_Phi_z()  const {
    Real z = zmax_ * math::cos(phases_.psi_z); // Kepler parametrize z(\psi_z)
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
void ghz::KerrBoundOrbit::compute_frequencies() {
    Real p3 = r3_*(one-e_)/M_;
    Real p4 = r4_*(one+e_)/M_;
    freqs_.f_t = get_T_r() + get_T_z() + a_*Lz_;  // Mino f_t = dt/d\lambda = Sigma dt/d\tau
    freqs_.f_phi = get_Phi_r() + get_Phi_z() - a_*Lz_;  // Mino f_phi = dphi/d\lambda = Sigma dphi/d\tau

    Real term1 = p_ - p3 - e_*(p_+p3*math::cos(phases_.psi_r));
    Real term2 = p_ - p4 + e_*(p_-p4*math::cos(phases_.psi_r));
    freqs_.f_r = M_ * math::sqrt(gamma_*term1*term2)/(one-e_*e_);  // (216) Pound and Wardell

    freqs_.f_z = a_*math::sqrt(
            gamma_ *
            (z1_*z1_-zmax_*zmax_*math::sqr(math::cos(phases_.psi_z)))
    );
}
/**
 * update_q_angles
 * @brief Update the mean action angles q_alpha
 *
 * - updates the private internal action angles q_alpha = Upsilon_alpha * lambda + q_alpha^0
 */
void ghz::KerrBoundOrbit::update_q_angles(ghz::Real mino_time_param) {
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
void ghz::KerrBoundOrbit::sample_T_and_Phi_for_fft(size_t Nr, size_t Nz)
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

void ghz::KerrBoundOrbit::sample_frequencies_for_fft(size_t Nr, size_t Nz) {

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

void ghz::KerrBoundOrbit::compute_q_from_psi(const std::vector<Real>& f, std::vector<Real>& q) {
    size_t N = f.size();

    std::vector<Real> invf(N), cum(N+1,0.0L);
    for (size_t i=0;i<N;i++) invf[i] = 1.0L / f[i];

    // cumulative integral
    for (size_t i=0;i<N;i++){
        Real dpsi = twoPi / N; // because ψ-grid is uniform
        cum[i+1] = cum[i] + 0.5L*(invf[i] + invf[(i+1)%N])*dpsi; // crude trapezoidal rule
    }
    // Normalize to get q = Upsilon * cum, rescaled to [0,2π)
    Real scale = twoPi / cum[N]; // = Upsilon / <f> automatically
    q.resize(N);
    for (size_t i=0;i<N;i++) q[i] = scale * cum[i]; // set q
}

void ghz::KerrBoundOrbit::compute_q_grids_from_samples(size_t Nr, size_t Nz,
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

    // Now compute q_r(ψ_r) and q_z(ψ_z) using your function:
    compute_q_from_psi(f_r_slice, q_r);
    compute_q_from_psi(f_z_slice, q_z);
}

// map FFT index -> integer mode k in range -(N/2)..(N/2-1)
inline int fft_index_to_mode(int idx, int N) {
    if (idx <= N/2) return idx;
    return idx - N;
}

void ghz::KerrBoundOrbit::compute_Delta_and_freq_modes(
        const std::vector<Complex>& samples_in,  // f(ψ) samples
        const std::vector<Real>& q_grid,         // q[ψ] sampled at same ψ-grid points <-- pass q_r_vals OR q_z_vals here
        Real Omega,                              // \Upsilon_α or Ω_α (torus freq)
        std::vector<Real>& Delta_out,            // Δψ or Δt or Δφ
        std::vector<Complex>& modes_out)         // Fourier modes c_k of frequency
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

    size_t N = samples_in.size();
    const Complex I(0.0L, 1.0L);

    Delta_out.assign(N, Real(0));
    modes_out.assign(N, Complex(0,0));

    // mean subtraction
    Complex mean = std::accumulate(samples_in.begin(), samples_in.end(),
                                   Complex(0.0L,0.0L)) / Real(N);

    std::vector<Complex> osc(N);
    for(size_t i=0;i<N;i++) osc[i] = samples_in[i] - mean;

    std::vector<Complex> raw_modes(N);

    fftwl_plan plan = fftwl_plan_dft_1d(int(N),
                                        reinterpret_cast<fftwl_complex*>(osc.data()),
                                        reinterpret_cast<fftwl_complex*>(raw_modes.data()),
                                        FFTW_FORWARD, FFTW_ESTIMATE);
    fftwl_execute(plan);
    fftwl_destroy_plan(plan);

    // normalize
    for(size_t k=0;k<N;k++) modes_out[k] = raw_modes[k] / Real(N);

    // construct Δ using *actual* q-grid
    int kmax = int(N/2);
    for(size_t i=0;i<N;i++){
        Real q = q_grid[i];
        Complex sum(0.0L,0.0L);

        for(int k=1; k<kmax; ++k){
            Complex ck = modes_out[k];
            Complex kC = Complex((Real)k,0.0L);

            sum += ck * std::exp(-I*kC*q) / (-I*kC*Omega);
            sum += std::conj(ck) * std::exp( I*kC*q) / ( I*kC*Omega);
        }

        Delta_out[i] = sum.real();
    }
}

// compute_Delta_psi_rz: build slices, call helper for f_r,f_z and for t,phi components
void ghz::KerrBoundOrbit::compute_Deltas(size_t Nr, size_t Nz)
{
    const Real two_pi = teuk::twoPi;
    // allocate slice containers with correct sizes
    std::vector<Complex> slice_r(Nr), slice_t_r(Nr), slice_phi_r(Nr);
    std::vector<Complex> slice_z(Nz), slice_t_z(Nz), slice_phi_z(Nz);

    // radial slice (choose j=0)
    for (size_t i = 0; i < Nr; ++i) {
        size_t idx = i * Nz + 0; // j=0 slice
        slice_r[i]     = Complex(f_samples_[idx].f_r, 0.0L);
        slice_t_r[i]   = Complex(f_samples_[idx].f_t, 0.0L);
        slice_phi_r[i] = Complex(f_samples_[idx].f_phi, 0.0L);
    }

    // polar slice (choose i=0)
    for (size_t j = 0; j < Nz; ++j) {
        size_t idx = 0 * Nz + j; // /i=0 slice
        slice_z[j]     = Complex(f_samples_[idx].f_z, 0.0L);
        slice_t_z[j]   = Complex(f_samples_[idx].f_t, 0.0L);
        slice_phi_z[j] = Complex(f_samples_[idx].f_phi, 0.0L);
    }

    // Prepare outputs (members public)
    Delta_psi_r_.resize(Nr); Delta_psi_z_.resize(Nz);
    Delta_t_r_.resize(Nr);   Delta_t_z_.resize(Nz);
    Delta_phi_r_.resize(Nr); Delta_phi_z_.resize(Nz);



    // compute radial deltas and modes
    compute_Delta_and_freq_modes(slice_r, qr_vals,  torus_freqs_.Ups_r, Delta_psi_r_, f_modes_.f_r_modes);
    compute_Delta_and_freq_modes(slice_t_r, qr_vals, torus_freqs_.Ups_r, Delta_t_r_,   f_modes_.t_r_modes);
    compute_Delta_and_freq_modes(slice_phi_r,qr_vals, torus_freqs_.Omega_r, Delta_phi_r_, f_modes_.phi_r_modes);

    // compute polar deltas and modes
    compute_Delta_and_freq_modes(slice_z, qz_vals,  torus_freqs_.Omega_z, Delta_psi_z_, f_modes_.f_z_modes);
    compute_Delta_and_freq_modes(slice_t_z, qz_vals, torus_freqs_.Omega_z, Delta_t_z_,   f_modes_.t_z_modes);
    compute_Delta_and_freq_modes(slice_phi_z, qz_vals, torus_freqs_.Omega_z, Delta_phi_z_, f_modes_.phi_z_modes);
}



void ghz::KerrBoundOrbit::compute_Delta_psi_rz(size_t Nr, size_t Nz)
{
    // FFT input arrays (real -> complex, you can also use complex input directly)
    // allocate 1D arrays for slices
    std::vector<Complex> f_r_samples(Nr);
    std::vector<Complex> f_z_samples(Nz);
    std::vector<Complex> f_r_modes(Nz);
    std::vector<Complex> f_z_modes(Nz);
    // Radial slice (i=0..Nr-1, pick j=0 slice for separable function using flattened index )
    for (size_t i = 0; i < Nr; ++i) {
        f_r_samples[i] = Complex(f_samples_[i*Nz].f_r, 0.0);
    }

// Subtract mean (oscillatory part)
    Real f_r_mean = 0.0;
    for (auto &v : f_r_samples) f_r_mean += v.real();
    f_r_mean /= Nr;
    for (auto &v : f_r_samples) v -= f_r_mean;

// Polar slice
    for (size_t j = 0; j < Nz; ++j) {
        f_z_samples[j] = Complex(f_samples_[j].f_z, 0.0);
    }

    // Subtract mean (oscillatory part)
    Real f_z_mean = 0.0;
    for (auto &v : f_z_samples) f_z_mean += v.real(); // add up all the samples
    f_z_mean /= Nz; // compute avg.
    for (auto &v : f_z_samples) v -= f_z_mean; // subtract mean from each sample

    fftwl_plan plan_r = fftwl_plan_dft_1d(Nr,
                                          reinterpret_cast<fftwl_complex*>(f_r_samples.data()),
                                          reinterpret_cast<fftwl_complex*>(f_r_modes.data()),
                                          FFTW_FORWARD, FFTW_ESTIMATE);

    fftwl_plan plan_z = fftwl_plan_dft_1d(Nz,
                                          reinterpret_cast<fftwl_complex*>(f_z_samples.data()),
                                          reinterpret_cast<fftwl_complex*>(f_z_modes.data()),
                                          FFTW_FORWARD, FFTW_ESTIMATE);

// execute plan
    fftwl_execute(plan_r);
    fftwl_execute(plan_z);
// Destroy plans
    fftwl_destroy_plan(plan_r);
    fftwl_destroy_plan(plan_z);

    std::vector<Real> delta_psi_r(Nr);
    std::vector<Real> delta_psi_z(Nz);

// Radial
    for (size_t i = 0; i < Nr; ++i) {
        Real q_r = qr_vals[i]; // twoPi * i / Nr;
        Complex sum = 0.0;
        for (int k = 1; k < Nr/2; ++k) {   // exclude k=0
            Complex kC = Complex(Real(k),0.0);
            sum += f_r_modes[k] * std::exp(-I*kC*q_r) / (-I*kC*torus_freqs_.Omega_r);
            sum += std::conj(f_r_modes[k]) * std::exp(I*kC*q_r) / (I*kC*torus_freqs_.Omega_r);
        }
        delta_psi_r[i] = sum.real();
    }

// Polar
    for (size_t j = 0; j < Nz; ++j) {
        Real q_z = qz_vals[j]; // twoPi * j / Nz;
        Complex sum = 0.0;
        for (int k = 1; k < Nz/2; ++k) {
            Complex kC = Complex(Real(k),0.0);
            sum += f_z_modes[k] * std::exp(-I*kC*q_z) / (-I*kC*torus_freqs_.Omega_z);
            sum += std::conj(f_z_modes[k]) * std::exp(I*kC*q_z) / (I*kC*torus_freqs_.Omega_z);
        }
        delta_psi_z[j] = sum.real();
    }
}