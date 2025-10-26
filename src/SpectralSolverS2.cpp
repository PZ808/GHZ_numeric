//
// Created by Peter Zimmerman on 26.10.25.
//

#include "../include/SpectralSolverS2.hpp"
#include <cmath>
#include <cassert>
#include "../include/GHPScalars.hpp"

namespace SpecS2 {

// ===========================================================
// Constructor: build nodes, weights, differentiation matrix
// ===========================================================
    SpectralSolver::SpectralSolver(int Nz, int Nphi)
            : Nz_(Nz), Nphi_(Nphi)
    {
        z_ = legendre_gauss_lobatto(Nz_);
        w_ = barycentric_weights();
        D_ = legendre_diff_matrix(z_);
    }

// ===========================================================
// dz derivative
// ===========================================================
    void SpectralSolver::dz(const vector<vector<Complex>>& f,
                            vector<vector<Complex>>& df) const
    {
        int Nphi = Nphi_;
        int Nz = Nz_;
        df.resize(Nz, vector<Complex>(Nphi));

        for(int j=0; j<Nphi; ++j) {
            Eigen::VectorXcd col(Nz), dcol(Nz);
            for(int i=0; i<Nz; ++i) col(i) = f[i][j];
            dcol = D_.cast<Complex>() * col;
            for(int i=0; i<Nz; ++i) df[i][j] = dcol(i);
        }
    }
    void SpectralSolver::dphi_fft(const std::vector<std::vector<Complex>>& f,
                                  std::vector<std::vector<Complex>>& df_dphi) const
    {
        int Nz   = Nz_;
        int Nphi = Nphi_;

        // Resize output if needed
        df_dphi.resize(Nz);
        for(int iz = 0; iz < Nz; ++iz)
            df_dphi[iz].resize(Nphi);

        std::vector<Complex> tmp_in(Nphi);   // temporary input copy for FFT
        std::vector<Complex> tmp_out(Nphi);  // FFT output

        for(int iz = 0; iz < Nz; ++iz)
        {
            // --- Step 0: copy the current φ-line ---
            tmp_in = f[iz];

            // --- Step 1: Forward FFT ---
            fftw_plan p_forward = fftw_plan_dft_1d(
                    Nphi,
                    reinterpret_cast<fftw_complex*>(tmp_in.data()),
                    reinterpret_cast<fftw_complex*>(tmp_out.data()),
                    FFTW_FORWARD, FFTW_ESTIMATE);

            fftw_execute(p_forward);
            fftw_destroy_plan(p_forward);

            // --- Step 2: Multiply by i*m for derivative in Fourier space ---
            for(int m = 0; m < Nphi; ++m)
            {
                int ms = (m <= Nphi / 2) ? m : m - Nphi; // handle negative frequencies
                tmp_out[m] *= Complex(0.0, static_cast<double>(ms));
            }

            // --- Step 3: Inverse FFT to φ-space ---
            fftw_plan p_backward = fftw_plan_dft_1d(
                    Nphi,
                    reinterpret_cast<fftw_complex*>(tmp_out.data()),
                    reinterpret_cast<fftw_complex*>(df_dphi[iz].data()),
                    FFTW_BACKWARD, FFTW_ESTIMATE);

            fftw_execute(p_backward);
            fftw_destroy_plan(p_backward);

            // --- Step 4: Normalize ---
            for(int k = 0; k < Nphi; ++k)
                df_dphi[iz][k] /= static_cast<double>(Nphi);
        }
    }
// ===========================================================
// Barycentric weights
// ===========================================================
    vector<double> SpectralSolver::barycentric_weights() const {
        int N = z_.size();
        vector<double> w(N,1.0);
        for(int j=0;j<N;++j){
            for(int k=0;k<N;++k) if(k!=j) w[j] /= (z_[j]-z_[k]);
        }
        return w;
    }

// ===========================================================
// Barycentric interpolation & derivative
// ===========================================================
    std::pair<Complex,Complex> SpectralSolver::barycentric_interp_and_derivative(
            const vector<Complex>& f, const vector<double>& w, double z0) const
    {
        int N = z_.size();
        Complex num(0,0), num_der(0,0);
        double den=0.0;

        for(int i=0;i<N;++i){
            double dz = z0 - z_[i];
            if(std::abs(dz)<1e-14) return {f[i], Complex(NAN,NAN)};
            double wi_dz = w[i]/dz;
            num += wi_dz*f[i];
            num_der += wi_dz*f[i]/dz;
            den += wi_dz;
        }

        Complex f_interp = num/den;
        Complex df_interp = num_der - f_interp*(num_der/num);
        return {f_interp, df_interp};
    }

// ===========================================================
// Legendre polynomial & derivative
// ===========================================================
    std::pair<double,double> SpectralSolver::legendre_P_and_dP(int n, double x){
        if(n==0) return {1.0,0.0};
        if(n==1) return {x,1.0};
        double Pnm1=1.0, Pn=x;
        for(int k=2;k<=n;++k){
            double Pnp1 = ((2*k-1)*x*Pn-(k-1)*Pnm1)/k;
            Pnm1 = Pn; Pn = Pnp1;
        }
        double dPn = n/(x*x-1)*(x*Pn-Pnm1);
        return {Pn,dPn};
    }

// ===========================================================
// Legendre Gauss-Lobatto nodes
// ===========================================================
    vector<double> SpectralSolver::legendre_gauss_lobatto(int N){
        vector<double> x(N);
        x[0] = -1.0; x[N-1] = 1.0;

        for(int i=1;i<N-1;++i){
            double xi=-cos(M_PI*i/(N-1));
            for(int it=0;it<20;++it){
                auto [P,dP] = legendre_P_and_dP(N-1,xi);
                double d2P = (2*xi*dP-(N-1)*N*P)/(1-xi*xi);
                double dx=-dP/d2P;
                xi+=dx;
                if(std::abs(dx)<1e-15) break;
            }
            x[i]=xi;
        }
        return x;
    }

// ===========================================================
// Legendre differentiation matrix
// ===========================================================
    MatrixXd SpectralSolver::legendre_diff_matrix(const vector<double>& x){
        int N=x.size();
        MatrixXd D=MatrixXd::Zero(N,N);
        vector<double> Pnm1(N);
        for(int i=0;i<N;++i) Pnm1[i]=legendre_P_and_dP(N-1,x[i]).first;

        for(int i=0;i<N;++i){
            for(int j=0;j<N;++j){
                if(i!=j) D(i,j)=Pnm1[i]/(Pnm1[j]*(x[i]-x[j]));
            }
        }
        for(int i=0;i<N;++i){
            double s=0.0;
            for(int j=0;j<N;++j) if(i!=j) s+=D(i,j);
            D(i,i)=-s;
        }
        return D;
    }

// ============================================================
// eth operator for a GHP scalar (spin-weighted derivative)
// ===========================================================
// @brief Computes the eth operator (spin-raising derivative) for a GHP scalar.
//
// Applies the spin-weighted derivative along the m-direction of the Newman-Penrose
// tetrad. The returned GHP scalar has its spin weight raised by one:
// (p,q) -> (p+1, q-1).
//
// @param f_in Input GHPScalar on the spectral grid (values stored in 2D vector).
// @return GHPScalar The eth derivative of the input, with updated spin weights.
//
    GHPField SpectralSolver::edth(const GHPField& f_in) const
    {
        int Nz = Nz_;
        int Nphi = Nphi_;

        // Allocate derivative arrays
        std::vector<std::vector<Complex>> df_dz(Nz, std::vector<Complex>(Nphi));
        std::vector<std::vector<Complex>> df_dphi(Nz, std::vector<Complex>(Nphi));

        // Copy input values
        const auto& f = f_in.values();

        // Compute spectral derivatives
        dz(f, df_dz);         // ∂/∂z
        dphi_fft(f, df_dphi); // ∂/∂φ

        // Compute eth combination
        std::vector<std::vector<Complex>> df_eth(Nz, std::vector<Complex>(Nphi));
        for(int i = 0; i < Nz; ++i) {
            double z = nodes()[i];
            double factor = std::sqrt(1.0 - z*z);
            for(int j = 0; j < Nphi; ++j) {
                df_eth[i][j] = -factor * (df_dz[i][j] + Complex(0.0,1.0) * df_dphi[i][j] / factor);
            }
        }

        // Return new field with raised GHP weights
        // Allocate new GHPField of same size
        GHPField out(Nz, Nphi, 0.0, f_in.p() + 1, f_in.q() - 1);

        // Fill in the computed values
        out.set_values(df_eth);
        return out;
    }

} // namespace SpecS2
