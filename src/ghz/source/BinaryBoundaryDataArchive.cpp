//
// Created by Peter Zimmerman on 20.03.26.
//
#include "ghz/source/BinaryBoundaryDataArchive.hpp"

#include <array>
#include <cstdint>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace ghz::source {

    namespace {
        constexpr std::array<char, 8> kMagic{{'G','H','Z','B','N','D','1','\0'}};
        constexpr std::uint32_t kVersion = 1u;

        template <typename T>
        T read_pod(std::ifstream& in, const char* what) {
            T x{};
            in.read(reinterpret_cast<char*>(&x), static_cast<std::streamsize>(sizeof(T)));
            if (!in) {
                throw std::runtime_error(std::string("BinaryBoundaryDataArchive: failed reading ") + what);
            }
            return x;
        }

        template <typename T>
        std::vector<T> read_pod_array(std::ifstream& in, std::size_t n, const char* what) {
            std::vector<T> out(n);
            if (n > 0) {
                in.read(reinterpret_cast<char*>(out.data()),
                        static_cast<std::streamsize>(n * sizeof(T)));
                if (!in) {
                    throw std::runtime_error(std::string("BinaryBoundaryDataArchive: failed reading array ") + what);
                }
            }
            return out;
        }

        bool axis_equal(const Axis& a, const Axis& b) {
            if (a.values.size() != b.values.size()) return false;
            for (std::size_t i = 0; i < a.values.size(); ++i) {
                if (!(a.values[i] == b.values[i])) return false;
            }
            return true;
        }

        bool same_metadata(const BoundaryMetadata& a, const BoundaryMetadata& b) {
            return (a.r_boundary == b.r_boundary) &&
                   (a.rp == b.rp) &&
                   (a.B  == b.B)  &&
                   (a.sym.negative_m_is_conjugate == b.sym.negative_m_is_conjugate) &&
                   (a.sym.y_storage == b.sym.y_storage);
        }
    } // namespace

    BinaryBoundaryDataArchive::BinaryBoundaryDataArchive(std::filesystem::path root_dir,
                                                         std::string prefix,
                                                         int probe_m,
                                                         PatchSide probe_patch)
            : root_dir_(std::move(root_dir)),
              prefix_(std::move(prefix))
    {
        LoadedBoundaryComponent probe = read_component_file_(BoundaryComponent::delta,
                                                             probe_m,
                                                             probe_patch);
        metadata_ = probe.meta;
    }

    BoundaryModeData BinaryBoundaryDataArchive::load_mode_patch(int m,
                                                                PatchSide patch) const
    {
        if (m >= 0) {
            return load_nonnegative_mode_patch_raw_(m, patch);
        }

        if (!metadata_.sym.negative_m_is_conjugate) {
            throw std::runtime_error(
                    "BinaryBoundaryDataArchive::load_mode_patch: negative m requested "
                    "but archive does not declare conjugation symmetry");
        }

        BoundaryModeData out = load_nonnegative_mode_patch_raw_(-m, patch);
        out.m = m;

        for (std::uint32_t craw = 0; craw < static_cast<std::uint32_t>(BoundaryComponent::count); ++craw) {
            auto& vec = out.components[craw];
            for (auto& v : vec) {
                v = std::conj(v);
            }
            if (patch != PatchSide::Left) {
                throw std::runtime_error(
                        "BinaryBoundaryDataArchive::load_mode_patch: only left boundary data is supported");
            }
        }

        return out;
    }

    BoundaryModeData BinaryBoundaryDataArchive::load_nonnegative_mode_patch_raw_(int m,
                                                                                 PatchSide patch) const
    {
        if (m < 0) {
            throw std::runtime_error(
                    "BinaryBoundaryDataArchive::load_nonnegative_mode_patch_raw_: m must be >= 0");
        }

        BoundaryModeData out;
        out.m = m;
        out.patch = patch;

        bool first = true;
        for (std::uint32_t craw = 0; craw < static_cast<std::uint32_t>(BoundaryComponent::count); ++craw) {
            const BoundaryComponent c = static_cast<BoundaryComponent>(craw);
            LoadedBoundaryComponent lc = read_component_file_(c, m, patch);

            ensure_same_metadata_(metadata_, lc.meta,
                                  "BinaryBoundaryDataArchive::load_nonnegative_mode_patch_raw_");

            if (first) {
                out.Y = lc.Y;
                first = false;
            } else {
                ensure_same_axes_(out.Y, lc.Y,
                                  "BinaryBoundaryDataArchive::load_nonnegative_mode_patch_raw_: Y mismatch");
            }

            out.component(c) = std::move(lc.data);
        }

        return out;
    }

    std::filesystem::path BinaryBoundaryDataArchive::component_path_(BoundaryComponent c,
                                                                     int m,
                                                                     PatchSide /*patch*/) const
    {
        std::ostringstream os;
        os << prefix_ << "_" << component_tag_(c)
           << "_m" << m
           << ".bin";
        return root_dir_ / os.str();
    }

    BinaryBoundaryDataArchive::LoadedBoundaryComponent
    BinaryBoundaryDataArchive::read_component_file_(BoundaryComponent c,
                                                    int m,
                                                    PatchSide patch) const
    {
        const std::filesystem::path path = component_path_(c, m, patch);

        std::ifstream in(path, std::ios::binary);
        if (!in) {
            throw std::runtime_error(
                    "BinaryBoundaryDataArchive: cannot open file " + path.string());
        }

        std::array<char, 8> magic{};
        in.read(magic.data(), static_cast<std::streamsize>(magic.size()));
        if (!in || magic != kMagic) {
            throw std::runtime_error(
                    "BinaryBoundaryDataArchive: bad magic in " + path.string());
        }

        const std::uint32_t version       = read_pod<std::uint32_t>(in, "version");
        const std::int32_t  m_file        = read_pod<std::int32_t>(in, "m");
        const std::uint32_t component_raw = read_pod<std::uint32_t>(in, "component");
        //const std::uint32_t patch_raw     = read_pod<std::uint32_t>(in, "patch");
        const std::uint32_t y_storage_raw = read_pod<std::uint32_t>(in, "y_storage");
        const std::uint32_t flags         = read_pod<std::uint32_t>(in, "flags");
        //const double        r_boundary_d  = read_pod<double>(in, "r_boundary");
        const double        rp_d          = read_pod<double>(in, "rp");
        const double        B_d           = read_pod<double>(in, "B");
        const std::uint64_t ny_u64        = read_pod<std::uint64_t>(in, "ny");

        if (version != kVersion) {
            throw std::runtime_error(
                    "BinaryBoundaryDataArchive: unsupported version in " + path.string());
        }

        if (component_raw >= static_cast<std::uint32_t>(BoundaryComponent::count)) {
            throw std::runtime_error(
                    "BinaryBoundaryDataArchive: bad component id in " + path.string());
        }

        const BoundaryComponent c_file = static_cast<BoundaryComponent>(component_raw);
        if (c_file != c) {
            throw std::runtime_error(
                    "BinaryBoundaryDataArchive: component mismatch in " + path.string());
        }

        if (m_file != m) {
            throw std::runtime_error(
                    "BinaryBoundaryDataArchive: m mismatch in " + path.string());
        }


        const PatchSide patch_file = PatchSide::Left;
        const std::size_t ny = static_cast<std::size_t>(ny_u64);
        if (ny == 0) {
            throw std::runtime_error(
                    "BinaryBoundaryDataArchive: empty axis in " + path.string());
        }

        const std::vector<double> Yd  = read_pod_array<double>(in, ny, "Y");
        const std::vector<double> buf = read_pod_array<double>(in, 2 * ny, "complex data");

        LoadedBoundaryComponent out;
        out.m = m_file;
        out.patch = patch_file;
        out.component = c_file;
        //out.meta.r_boundary = Real(r_boundary_d);
        out.meta.rp = Real(rp_d);
        out.meta.B  = Real(B_d);
        out.meta.sym = symmetry_from_fields_(flags, y_storage_raw);

        out.Y.values.resize(ny);
        for (std::size_t j = 0; j < ny; ++j) {
            out.Y.values[j] = Real(Yd[j]);
        }

        out.data.resize(ny);
        for (std::size_t k = 0; k < ny; ++k) {
            out.data[k] = Complex(Real(buf[2 * k]), Real(buf[2 * k + 1]));
        }

        return out;
    }
    std::string BinaryBoundaryDataArchive::component_tag_(BoundaryComponent c)
    {
        switch (c) {
            case BoundaryComponent::delta:      return "ll_delta";
            case BoundaryComponent::deltaPrime: return "ll_deltaPrime";
            default:
                throw std::runtime_error("BinaryBoundaryDataArchive: unknown component");
        }
    }

    void BinaryBoundaryDataArchive::ensure_same_metadata_(const BoundaryMetadata& a,
                                                          const BoundaryMetadata& b,
                                                          const std::string& where)
    {
        if (!same_metadata(a, b)) {
            throw std::runtime_error(where + ": metadata mismatch across files");
        }
    }

    void BinaryBoundaryDataArchive::ensure_same_axes_(const Axis& a,
                                                      const Axis& b,
                                                      const std::string& where)
    {
        if (!axis_equal(a, b)) {
            throw std::runtime_error(where);
        }
    }

    SourceSymmetry BinaryBoundaryDataArchive::symmetry_from_fields_(std::uint32_t flags,
                                                                    std::uint32_t y_storage_raw)
    {
        SourceSymmetry sym;
        sym.negative_m_is_conjugate = (flags & 0x1u) != 0u;

        switch (y_storage_raw) {
            case 0u:
                sym.y_storage = YStorage::Full;
                break;
            case 1u:
                sym.y_storage = YStorage::HalfWithEvenReflection;
                break;
            default:
                throw std::runtime_error("BinaryBoundaryDataArchive: unsupported y_storage id");
        }

        return sym;
    }

} // namespace ghz::source