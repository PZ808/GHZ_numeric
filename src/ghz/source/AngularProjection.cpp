//
// AngularProjection.cpp
//

#include "ghz/source/AngularProjection.hpp"

#include <sstream>
#include <stdexcept>

namespace ghz::source {

Complex project_angular_slice_z(const std::vector<Complex>& source_values,
                                const AngularModeSamples& harmonic,
                                const std::vector<Real>& weights_z) {
    if (source_values.size() != harmonic.values.size() ||
        source_values.size() != weights_z.size()) {
        std::ostringstream os;
        os << "project_angular_slice_z: size mismatch among source, harmonic, "
              "and weights";
        throw std::runtime_error(os.str());
    }

    Complex sum = teuk::zeroC;
    for (size_t iz = 0; iz < source_values.size(); ++iz) {
        sum += weights_z[iz] * source_values[iz] * std::conj(harmonic.values[iz]);
    }

    return teuk::twoPi * sum;
}

std::vector<Complex> project_angular_mode(const TeukolskyField& source,
                                          const AngularModeSamples& harmonic,
                                          const std::vector<Real>& weights_z) {
    if (source.Nz() != harmonic.values.size() ||
        source.Nz() != weights_z.size()) {
        std::ostringstream os;
        os << "project_angular_mode: angular sample count mismatch";
        throw std::runtime_error(os.str());
    }

    std::vector<Complex> out(source.Nr(), teuk::zeroC);
    std::vector<Complex> row(source.Nz(), teuk::zeroC);

    for (size_t ir = 0; ir < source.Nr(); ++ir) {
        for (size_t iz = 0; iz < source.Nz(); ++iz) {
            row[iz] = source(ir, iz).value();
        }
        out[ir] = project_angular_slice_z(row, harmonic, weights_z);
    }

    return out;
}

} // namespace ghz::source
