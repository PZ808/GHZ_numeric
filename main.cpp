#include "ghz/core/GhzTypes.hpp"
#include "ghz/geom/KerrParams.hpp"
#include "ghz/geom/KerrMetric.hpp"
#include "ghz/geom/KerrMetricOutgoing.hpp"

#include "ghz/spectral/SpectralDiffer.hpp"
#include "ghz/spectral/SpectralCoordinateMaps.hpp"

#include "ghz/geom/DataDomain.hpp"

#include "ghz/orbit/KerrOrbit.hpp"
#include "ghz/source/BinaryEffectiveSourceArchive.hpp"
#include "ghz/source/ConditionSource.hpp"
#include "ghz/source/BinaryBoundaryDataArchive.hpp"
#include "ghz/source/ConditionBoundaryData.hpp"

#include "ghz/transport/Corrector.hpp"
#include "ghz/transport/collocation/XmmbarWorldtubeSolve.hpp"

#include <filesystem>
#include <iostream>
#include <stdexcept>

int main()
{
    using teuk::Real;
    using teuk::Complex;
    using namespace ghz::transport;
    using namespace ghz::source;
    std::cerr << "start main" << std::endl;

    const Real M = 1.0;
    const Real a = 0.9;

    KerrParams params(M, a);
    KerrMetric kerr(params);

    // create an orbit object
    try {

        // Define parameters for KerrBoundOrbit
        Real p = 15.0;       // Example value for semi-latus rectum
        Real e = 0.5;        // Example value for eccentricity
        Real inc = 0.0; // Example value for inclination (45 degrees)
        signed int chi = 1;  // Prograde orbit
        size_t Nr = 101;     // Odd number of radial samples
        size_t Nz = 101;     // Odd number of polar samples

        // Initialize KerrBoundOrbit object
        orbit::KerrBoundOrbit equatorial_orbit(kerr, p, e, inc, chi, Nr, Nz);

        equatorial_orbit.initialize_orbit();

        // Export FFT data to a file
        equatorial_orbit.export_equatorial_fft_data("fft_data.csv");

        // Export FFT metadata to a file
        equatorial_orbit.export_equatorial_fft_metadata("fft_metadata.txt");

        std::cout << "FFT data and metadata exported successfully.\n";
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
    }


    try {
        // ------------------------------------------------------------
        // Basic Kerr setup
        // ------------------------------------------------------------

        KerrMetricOutgoing kerr_out(params, kerr);

        // ------------------------------------------------------------
        // Mode / field metadata
        // ------------------------------------------------------------
        const int m_mode  = 0;
        const int kr_mode = 0;
        const int kz_mode = 0;

        const Modes modes{m_mode, kr_mode, kz_mode};
        const GHPType out_type{0, 0};

        // ------------------------------------------------------------
        // Worldtube split
        // ------------------------------------------------------------
        const Real r_min = 8.0;
        const Real r_p   = 10.0;
        const Real r_max = 12.0;

        ghz::numeric::PunctureTwoDomainSplit split(r_min, r_p, r_max);

        // ------------------------------------------------------------
        // Spectral resolutions
        // ------------------------------------------------------------
        const std::size_t Nz       = 41;
        const std::size_t Nr_left  = 31;
        const std::size_t Nr_right = 31;

        ::spectral::SpectralDiffer differ_left (Nz, Nr_left);
        ::spectral::SpectralDiffer differ_right(Nz, Nr_right);

        const RVector z_grid = differ_left.lgl_nodes();

        // ------------------------------------------------------------
        // Archive locations
        // ------------------------------------------------------------
        const std::filesystem::path src_dir =
                "/Users/antares/Projects/physics_codes/GHZ_numeric/Data/TeffBinary";

        std::cerr << "constructed src archive" << std::endl;
        const std::filesystem::path bdy_dir =
                "/Users/antares/Projects/physics_codes/GHZ_numeric/Data/TeffBoundaryBinary";
        std::cerr << "constructed bdy archive" << std::endl;
        // ------------------------------------------------------------
        // Bulk source archive + conditioned sources
        // ------------------------------------------------------------
        BinaryEffectiveSourceArchive src_archive(
                src_dir, "Teff", m_mode, PatchSide::Left);

        ConditionSource src_left (src_archive, m_mode, PatchSide::Left);
        std::cerr << "constructed conditioned source left" << std::endl;
        ConditionSource src_right(src_archive, m_mode, PatchSide::Right);
        std::cerr << "constructed conditioned source right" << std::endl;

        // ------------------------------------------------------------
        // Boundary archive + conditioned boundary data
        // ------------------------------------------------------------
        BinaryBoundaryDataArchive bdy_archive( bdy_dir, "TeffBdy", m_mode, PatchSide::Left);
        std::cerr << "built boundary archive" << std::endl;

        ConditionBoundaryData bdy_left(bdy_archive, m_mode, PatchSide::Left);
        std::cerr << "constructed conditioned boundary data" << std::endl;

        // ------------------------------------------------------------
        // Solve X_{m \bar m} on the worldtube
        // ------------------------------------------------------------
        std::cerr << "calling solve_xmmbar_worldtube" << std::endl;
        TwoDomainField xmmbar = solve_xmmbar_worldtube(
                differ_left,
                differ_right,
                split,
                kerr_out,
                z_grid,
                modes,
                out_type,
                src_left,
                src_right,
                bdy_left
        );

        // ------------------------------------------------------------
        // Quick diagnostics
        // ------------------------------------------------------------
        std::cout << "Solved X_mmbar worldtube.\n";
        std::cout << "left field:  Nr=" << xmmbar.left.Nr()
                  << " Nz=" << xmmbar.left.Nz() << "\n";
        std::cout << "right field: Nr=" << xmmbar.right.Nr()
                  << " Nz=" << xmmbar.right.Nz() << "\n";

        const std::size_t iz_mid = z_grid.size() / 2;
        const std::size_t iL0 = 0;
        const std::size_t iLp = xmmbar.left.Nr() - 1;
        const std::size_t iR0 = 0;
        const std::size_t iR1 = xmmbar.right.Nr() - 1;

        std::cout << "z_mid = " << z_grid[iz_mid] << "\n";
        std::cout << "Xmmbar.left(r_min, z_mid) = "
                  << xmmbar.left(iL0, iz_mid).value() << "\n";
        std::cout << "Xmmbar.left(r_p-, z_mid)  = "
                  << xmmbar.left(iLp, iz_mid).value() << "\n";
        std::cout << "Xmmbar.right(r_p+, z_mid) = "
                  << xmmbar.right(iR0, iz_mid).value() << "\n";
        std::cout << "Xmmbar.right(r_max, z_mid)= "
                  << xmmbar.right(iR1, iz_mid).value() << "\n";

        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "main failed: " << e.what() << "\n";
        return 1;
    }
}