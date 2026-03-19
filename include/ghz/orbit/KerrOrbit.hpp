//
// Created by Peter Zimmerman on 31.10.25.
//

#ifndef GHZ_NUMERIC_KERRORBIT_HPP
#define GHZ_NUMERIC_KERRORBIT_HPP

#include "ghz/core/GhzTypes.hpp"
#include "ghz/geom/KerrParams.hpp"
#include "ghz/core/MathMacros.hpp"
#include "ghz/geom/KerrMetric.hpp"
#include "ghz/geom/Coords.hpp"
#include "ghz/special/EllipticIntegrals.hpp"
#include "SourceTorusFrequencies.hpp"
#include "Splines.hpp"
#include <fftw3.h>
#include <cstddef>
#include <ranges>

namespace orbit {

    using Real = teuk::Real;
    using vectorR = std::vector<Real>;
    using vectorC = std::vector<Complex>;

    class KerrOrbitBase {
    protected:
        const KerrMetric& gKerr;
        Real a_, M_;

    public:
        explicit KerrOrbitBase(const KerrMetric& km) : gKerr(km)  {
            a_ = gKerr.a();
            M_ = gKerr.M();
        }

        virtual ~KerrOrbitBase() = default;

        [[nodiscard]] Real a() const { return a_; }
        [[nodiscard]] Real M() const { return M_; }

    };

/**
*  General Bound Kerr Orbit (action–angle representation)
*/
    class KerrBoundOrbit : public KerrOrbitBase {
    protected:

    protected:
        struct TorusFrequencies {
            Real Ups_t, Ups_phi, Ups_r, Ups_z; // avg frequencies dq_\mu/d\lambda
            Real Omega_t, Omega_phi, Omega_r, Omega_z; // avg frequencies dq_\mu/dt
        };
        struct Actions {
            Real J_t;    // azimuthal action (= Lz)
            Real J_phi;    // azimuthal action (= Lz)
            Real J_r;      // radial action
            Real J_z;      // polar action
        };
        // uniform advance angles dq/d\lambda = constant
        struct TorusAngles {
            Real q_t, q_phi;
            Real q_r, q_z;
        };
        // non-uniform advance angles d\psi/d\lambda = non constant f
        struct Phases {
            Real psi_t, psi_phi;
            Real psi_r, psi_z;
        };

        // instantaneous frequencies a = r,z
        struct Frequencies_fa {
            Real f_r;
            Real f_z;
        };
        // instantaneous frequencies
        struct Frequencies {
            Real f_t;
            Real f_phi;
            Real f_r;
            Real f_z;
        };
        // corresponding Fourier coefficients (after 1-d FFT)
        struct FrequencyModes {
            // f_a (k-modes)
            vectorC f_r_modes;
            vectorC f_z_modes;
            // f_t = T_r + T_z + a L_z (k-modes)
            vectorC T_r_modes;
            vectorC T_z_modes;
            // f_phi = Phi_r + Phi_z - a E (k-modes)
            vectorC Phi_r_modes;
            vectorC Phi_z_modes;
        };

        vectorR qr_vals;
        vectorR qz_vals;
        vectorR psi_r_vals;
        vectorR psi_z_vals;
        vectorR q_of_psi_r_vals;  // q_r(psi_r_i) at uniform psi_r_i
        vectorR q_of_psi_z_vals;  // q_z(psi_z_j) at uniform psi_z_j


    private:
        Real p_;        // semi-latus
        Real e_;        // eccentricity
        Real rp_, ra_;  // apoapsis and periaspsis where R(r) = 0
        Real r3_, r4_;  // remaining roots of R(r) = 0
        Real inc_;      // inclination of the orbit relative to the equatorial plane
        Real zmax_, zmin_;   // zmax = \sin(inc) where  \Theta(z) = 0
        Real z1_, z2_;       // unphysical roots z2 =  zmax
        signed int  chi_;   //  prograde or retrograde
        size_t Nz_, Nr_;

        std::vector<Frequencies> f_samples_;
        std::vector<Phases> psi_samples_;
        vectorR T_r_samples_, T_z_samples_, Phi_r_samples_, Phi_z_samples_;
        // modes containers (store normalized modes c_k = raw/N)
        vectorC psi_r_modes_, psi_z_modes_, t_r_modes_, t_z_modes_, phi_r_modes_, phi_z_modes_;
        vectorR Delta_psi_r_, Delta_psi_z_, Delta_t_r_, Delta_t_z_, Delta_phi_r_, Delta_phi_z_;

        Actions actions_;
        TorusFrequencies mino_torus_freqs, torus_freqs_; // averaged frequencies
        Frequencies freqs_; // instantaneous frequencies
        TorusAngles torus_angles_;     // instantaneous angles
        Phases phases_;     // instantaneous phases
        FrequencyModes f_modes_;

        Real E2_, E_, Lz_, Q_;
        Real alpha_, beta_, gamma_;

        // Determinant system builder using d_, f_, g_, h_
        struct Determinants_ {
            Real dg, gd;
            Real dh, hd;
            Real gh, hg;
            Real fh, hf;
            Real fg, gf;
        };

        // helpers
        inline const Real r3() const {
            return half*( alpha_+sqrt(math::sqr(alpha_)-Real(4.0)*beta_) ) ; }
        inline const Real r4() const { return beta_/r3(); }

        inline const Real d_(Real r) const {
            return (gKerr.Delta(r)*(r*r + gKerr.a()*gKerr.a()*zmax_)) /
                   (gKerr.M()*gKerr.M()* gKerr.M()*gKerr.M());
        }
        inline const Real f_(Real r) const {
            return ( r*r*r*r + gKerr.a()*gKerr.a() * (r*(r+Real(2.0)) + zmax_*zmax_*gKerr.Delta(r)) ) /
                   (gKerr.M()*gKerr.M()* gKerr.M()*gKerr.M());
        }
        inline const Real g_(Real r) const { return (Real(2.0)*gKerr.a()*r) / (gKerr.M()*gKerr.M());
        }
        inline const Real h_(Real r) const {
            return ( r*(r-Real(2.0)*gKerr.M()) + zmax_*zmax_*gKerr.Delta(r)/(1-zmax_*zmax_) ) /
                   (gKerr.M()*gKerr.M()) ;
        }
        // helper functions appearing the mino frequencies \Upsilon_\mu
        inline const Real G_(const Real& ka) const {
            return EllipticIntegrals::computeSecondKind(ka)/EllipticIntegrals::computeFirstKind(ka);
        }
        inline const Real F_pm_(const Real & r_plus_or_minus, const Real& kr) const {
            Real h = (ra_-rp_)*(r3_-r_plus_or_minus)/((ra_-r3_)*(rp_-r_plus_or_minus));
            return (rp_-r3_)*EllipticIntegrals::computeThirdKind(h,kr)
                   /EllipticIntegrals::computeFirstKind(kr);
        }
        inline const Real Fr_(const Real& kr) const {
            Real h =(ra_-rp_)/(ra_-r3_);
            return (rp_-r3_)*EllipticIntegrals::computeThirdKind(h,kr)
                   /EllipticIntegrals::computeFirstKind(kr);
        }
        inline Determinants_ computeDeterminants_(const Real& rp, const Real& ra) const {
            Determinants_ det{};

            // Evaluate functions at both radii
            Real dp = d_(rp), fp = f_(rp), gp = g_(rp), hp = h_(rp);
            Real da = d_(ra), fa = f_(ra), ga = g_(ra), ha = h_(ra);
            // Determinants for system solving (example structure)

            auto det2 = [](const Real& x1, const Real& x2,
                           const Real& y1, const Real& y2) { return x1 * y2 - x2 * y1; };
            det.dg = det2(dp, da, gp, ga);
            det.gd = -det.dg;
            det.dh = det2(dp, da, hp, ha);
            det.hd = -det.dh;
            det.gh = det2(gp, ga, hp, ha);
            det.hg = -det.gh;
            det.fh = det2(fp, fa, hp, ha);
            det.hf = -det.fh;
            det.fg = det2(fp, fa, gp, ga);
            det.gf = -det.fg;
            return det;
        } //computeDeterminants

        // --- helper: build all splines ---
        void build_splines() {
            interp_dpsi_r  = PeriodicSpline(qr_vals, Delta_psi_r_);
            interp_phi_r   = PeriodicSpline(qr_vals, Delta_phi_r_);
            interp_t_r     = PeriodicSpline(qr_vals, Delta_t_r_);

            interp_dpsi_z  = PeriodicSpline(qz_vals, Delta_psi_z_);
            interp_phi_z   = PeriodicSpline(qz_vals, Delta_phi_z_);
            interp_t_z     = PeriodicSpline(qz_vals, Delta_t_z_);
        }
        std::vector<Complex> resample_to_uniform_q(
                const std::vector<Complex>& samples_in,
                const std::vector<Real>& q_grid) const;

    protected:
        fftwl_plan fft_plan_r_{nullptr}; // for Nr
        fftwl_plan fft_plan_z_{nullptr}; // for Nz
        fftwl_complex* fft_in_r_{nullptr};
        fftwl_complex* fft_out_r_{nullptr};
        fftwl_complex* fft_in_z_{nullptr};
        fftwl_complex* fft_out_z_{nullptr};

        size_t fft_size_r_{};
        size_t fft_size_z_{};

    public:

        // Constructor
        /*
         *  The workflow is:
         * 1. Sample f_alpha(\psi_r,\psi_z) on a uniform grid of (psi_r,psi_z)
         *   phases. This is what sample_frequencies_for_fft() does
         * 2. Use FFT to extract the 1D fourier coefficients f_z^{k}, f_r^{k}, T_r^_{k}, T_z^{k}, Phi_r^{k}, Phi_z^{k}
         * where T_r(r) is the separated part of the f_t frequency depending on r etc,
         * 3. Use (251) to (253) to compute the oscillating pieces
         * Delta psi_alpha using the FFTs and mean frequencies Upsilon
         * 4. Construct the oscillating part of phases and the mean angles q(psi)
         */
        KerrBoundOrbit(const KerrMetric& gKerr,
                       Real p, Real e, Real inc, signed int chi,
                       size_t Nr, size_t Nz)
                : KerrOrbitBase(gKerr),
                  p_(p), e_(e), inc_(inc), chi_(chi),
                  Nr_(Nr), Nz_(Nz),
                  psi_r_vals(Nr), psi_z_vals(Nz)
        {
            assert((e_ >= 0.0L && e_ < 1.0L) &&
                   "Error: Eccentricity must be in [0,1) for bound orbits!");
            assert((inc_ >= 0.0L && inc_ <= M_PI) &&
                   "Error: Inclination must be in [0, pi]!");
            assert((std::abs(chi_) == 1) &&
                   "Error: chi must be +1 (prograde) or -1 (retrograde)!");
            assert((Nr_ % 2 == 1) && "Odd Nr is preferred to avoid Nyquist ambiguity.");
            assert((Nz_ % 2 == 1) && "Odd Nz is preferred to avoid Nyquist ambiguity.");

            // physical roots and turning points
            // compute turning points in Keplerian parametrization
            zmax_  = abs(sin(inc_));
            zmin_ = -zmax_;
            z2_ =  zmax_; //  polar turning point
            rp_ = p_ * M_ / (Real(1.0) + e_);  // periapsis radial turning point
            ra_  = p_ * M_/ (Real(1.0) - e_); // apoapsis radial turning point

            // set constants of motion
            set_constants_of_motion(); // computes E_, Lz_, Q_, alpha_, beta_, gamma_;
            // roots (unphysical, but needed for freq computations)
            z1_ = sqrt(Q_/(a_*a_*gamma_*zmax_*zmax_));

            r3_ = r3();
            r4_ = r4();
            assert(r4_ < r3_ && r3_ < rp_ && rp_ < ra_ &&
                   "Error: Root ordering violation!");
            assert(z2_ < z1_ &&
                   "Error: Polar root ordering violation!");

            // print out roots for debugging
            std::cout << "KerrBoundOrbit initialized with roots: \n"
                      << " rp = " << rp_ << ", ra = " << ra_
                      << ", r3 = " << r3_ << ", r4 = " << r4_ << "\n"
                      << " z1 = " << z1_ << ", z2 = " << z2_ << ", zmax = " << zmax_ << "\n";

            compute_torus_frequencies(); // compute Upsilons and Omega (avg frequencies in mino and BL time)

            // initialize to periapsis and max z
            torus_angles_ = {0.0, 0.0, 0.0, 0.0};
            phases_ = {0.0, 0.0, 0.0, 0.0};

            // initialize psi grids
            for (size_t i = 0; i < Nr; ++i) psi_r_vals[i] = 2.0*M_PI*i / Nr;
            for (size_t j = 0; j < Nz; ++j) psi_z_vals[j] = 2.0*M_PI*j / Nz;

            prepare_fft_plans();

        } // constructor KerrBoundOrbit

        ~KerrBoundOrbit() {
            if (fft_plan_r_) fftwl_destroy_plan(fft_plan_r_);
            if (fft_in_r_)  fftwl_free(fft_in_r_);
            if (fft_out_r_) fftwl_free(fft_out_r_);

            if (fft_plan_z_) fftwl_destroy_plan(fft_plan_z_);
            if (fft_in_z_)  fftwl_free(fft_in_z_);
            if (fft_out_z_) fftwl_free(fft_out_z_);
        }
        void initialize_orbit();
        void prepare_fft_plans();
        void free_fft();                     // cleanup fft
        // -------------------------------------------------
        // Computation methods (symbolic/numerical stubs)
        // -------------------------------------------------
        void set_constants_of_motion();
        void compute_torus_frequencies();
        void compute_frequencies();

        inline FrequencyModes get_freq_modes() { return f_modes_; }

        void update_q_angles(const Real& mino_time_param);
        BLCoords eval_at_Mino(const Real& lambda) const;

        void sample_frequencies_for_fft(size_t Nr, size_t Nz); // f_a
        void compute_q_grids_from_samples(const size_t& Nr, const size_t& Nz,
                                          std::vector<Real> &qr_vals, std::vector<Real> &qz_vals);
        void compute_Deltas(size_t Nr, size_t Nz);
        void compute_q_from_psi(const std::vector<Real>& f, std::vector<Real>& q);

        // -------------------------------------------------
        // Accessors
        // -------------------------------------------------
        [[nodiscard]] const Actions& actions() const { return actions_; }
        [[nodiscard]] const Frequencies& freqs() const { return freqs_; }
        [[nodiscard]] const Phases& angles() const { return phases_; }

        TorusFrequencies get_torus_frequencies() const { return torus_freqs_; } // return \Omega_t, \Omega_phi, \Omega_r, \Omega_z
        TorusFrequencies get_mino_torus_frequencies() const { return mino_torus_freqs; }

        Real get_T_r() const;
        [[nodiscard]] Real get_T_z() const;

        [[nodiscard]] Real get_Phi_r() const;
        [[nodiscard]] Real get_Phi_z() const;

        std::vector<Real> get_psi_r() const;
        std::vector<Real> get_psi_z() const;
        std::vector<Real> get_psi_t() const;
        std::vector<Real> get_psi_phi() const;

        // expose the torus source frequencies
        [[nodiscard]] orbit::SourceFrequencies source_frequencies() const noexcept {
            return { torus_freqs_.Omega_r, torus_freqs_.Omega_z, torus_freqs_.Omega_phi };
        }

        // --- splines ---
        PeriodicSpline interp_dpsi_r, interp_phi_r, interp_t_r;
        PeriodicSpline interp_dpsi_z, interp_phi_z, interp_t_z;



        void compute_Delta_and_freq_modes_SIMD(const std::vector<Complex> &samples_in,
                                          const std::vector<Real>& q_grid,
                                          const Real& Omega,
                                          std::vector<Real> &Delta_out,
                                          std::vector<Complex> &modes_out);


        void export_trajectory_stream(const std::string& filename, const Real& lambda_max,
                                      const Real& dlambda, size_t output_stride = 1);

    };  //  class KerrBoundOrbit

    class KerrCircularEquatorialOrbit : public KerrOrbitBase {
    private:
        Real r0_;
        int chi_;
        Real E_;
        Real Lz_;
        Real Ups_t_;
        Real Ups_phi_;
        Real Omega_phi_;

    public:
        KerrCircularEquatorialOrbit(const KerrMetric& km, Real r0, int chi);

        [[nodiscard]] Real energy() const noexcept { return E_; }
        [[nodiscard]] Real angular_momentum() const noexcept { return Lz_; }
        [[nodiscard]] Real upsilon_t() const noexcept { return Ups_t_; }
        [[nodiscard]] Real upsilon_phi() const noexcept { return Ups_phi_; }
        [[nodiscard]] Real Omega_phi() const noexcept { return Omega_phi_; }

        BLCoords eval_at_Mino(const Real& lambda) const {
            return BLCoords{
                    Ups_t_ * lambda,
                    r0_,
                    Real(0.0),   // z = cos(theta) = 0 on the equator
                    Ups_phi_ * lambda
            };
        }

        [[nodiscard]] orbit::SourceFrequencies source_frequencies() const noexcept {
            return { Real(0.0), Real(0.0), Omega_phi_ };
        }
    };

} // namespace orbit

#endif //GHZ_NUMERIC_KERRORBIT_HPP
