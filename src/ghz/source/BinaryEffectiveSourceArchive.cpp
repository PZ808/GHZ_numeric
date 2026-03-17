//
// Created by Peter Zimmerman on 17.03.26.
//
#include "ghz/source/BinaryEffectiveSourceArchive.hpp"

#include <array>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace ghz::source {

    namespace {
        constexpr std::array<char, 8> kMagic{{'G','H','Z','S','R','C','1','\0'}};
        constexpr std::uint32_t kVersion = 1u;

        template <typename T>
        T read_pod(std::ifstream& in, const char* what) {
            T x{};
            in.read(reinterpret_cast<char*>(&x), static_cast<std::streamsize>(sizeof(T)));
            if (!in) {
                throw std::runtime_error(std::string("BinaryEffectiveSourceArchive: failed reading ") + what);
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
                    throw std::runtime_error(std::string("BinaryEffectiveSourceArchive: failed reading array ") + what);
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

        bool same_metadata(const SourceMetadata& a, const SourceMetadata& b) {
            return (a.rp == b.rp) &&
                   (a.B  == b.B)  &&
                   (a.sym.negative_m_is_conjugate == b.sym.negative_m_is_conjugate) &&
                   (a.sym.y_storage == b.sym.y_storage);
        }
    } // namespace

    BinaryEffectiveSourceArchive::BinaryEffectiveSourceArchive(std::filesystem::path root_dir,
                                                               std::string prefix,
                                                               int probe_m,
                                                               PatchSide probe_patch)
            : root_dir_(std::move(root_dir)),
              prefix_(std::move(prefix))
    {
        LoadedComponent probe = read_component_file_(Component::ll, probe_m, probe_patch);
        metadata_ = probe.meta;
    }

    ModeData BinaryEffectiveSourceArchive::load_nonnegative_mode_patch_raw_(int m,
                                                                            PatchSide patch) const {
        if (m < 0) {
            throw std::runtime_error(
                    "BinaryEffectiveSourceArchive::load_nonnegative_mode_patch_raw_: m must be >= 0"
            );
        }

        ModeData out;
        out.m = m;
        out.patch = patch;

        bool first = true;
        for (std::uint32_t craw = 0; craw < static_cast<std::uint32_t>(Component::count); ++craw) {
            const Component c = static_cast<Component>(craw);
            LoadedComponent lc = read_component_file_(c, m, patch);

            ensure_same_metadata_(metadata_, lc.meta,
                                  "BinaryEffectiveSourceArchive::load_nonnegative_mode_patch_raw_");

            if (first) {
                out.X = lc.X;
                out.Y = lc.Y;
                first = false;
            } else {
                ensure_same_axes_(out.X, lc.X,
                                  "BinaryEffectiveSourceArchive::load_nonnegative_mode_patch_raw_: X mismatch");
                ensure_same_axes_(out.Y, lc.Y,
                                  "BinaryEffectiveSourceArchive::load_nonnegative_mode_patch_raw_: Y mismatch");
            }

            out.component(c) = std::move(lc.data);
        }

        return out;
    }

    std::filesystem::path BinaryEffectiveSourceArchive::component_path_(Component c,
                                                                        int m,
                                                                        PatchSide patch) const {
        std::ostringstream os;
        os << prefix_ << "_" << component_tag_(c)
           << "_m" << m
           << "_" << patch_tag(patch)
           << ".bin";
        return root_dir_ / os.str();
    }

    BinaryEffectiveSourceArchive::LoadedComponent
    BinaryEffectiveSourceArchive::read_component_file_(Component c,
                                                       int m,
                                                       PatchSide patch) const
    {
        const std::filesystem::path path = component_path_(c, m, patch);

        std::ifstream in(path, std::ios::binary);
        if (!in) {
            throw std::runtime_error(
                    "BinaryEffectiveSourceArchive: cannot open file " + path.string()
            );
        }

        // magic
        std::array<char, 8> magic{};
        in.read(magic.data(), static_cast<std::streamsize>(magic.size()));
        if (!in || magic != kMagic) {
            throw std::runtime_error(
                    "BinaryEffectiveSourceArchive: bad magic in " + path.string()
            );
        }

        const std::uint32_t version       = read_pod<std::uint32_t>(in, "version");
        const std::int32_t  m_file        = read_pod<std::int32_t>(in,  "m");
        const std::uint32_t component_raw = read_pod<std::uint32_t>(in, "component");
        const std::uint32_t patch_raw     = read_pod<std::uint32_t>(in, "patch");
        const std::uint32_t y_storage_raw = read_pod<std::uint32_t>(in, "y_storage");
        const std::uint32_t flags         = read_pod<std::uint32_t>(in, "flags");
        const double        rp_d          = read_pod<double>(in, "rp");
        const double        B_d           = read_pod<double>(in, "B");
        const std::uint64_t nx_u64        = read_pod<std::uint64_t>(in, "nx");
        const std::uint64_t ny_u64        = read_pod<std::uint64_t>(in, "ny");

        if (version != kVersion) {
            throw std::runtime_error(
                    "BinaryEffectiveSourceArchive: unsupported version in " + path.string()
            );
        }

        if (component_raw >= static_cast<std::uint32_t>(Component::count)) {
            throw std::runtime_error(
                    "BinaryEffectiveSourceArchive: bad component id in " + path.string()
            );
        }

        const Component c_file = static_cast<Component>(component_raw);
        if (c_file != c) {
            throw std::runtime_error(
                    "BinaryEffectiveSourceArchive: component mismatch in " + path.string()
            );
        }

        if (m_file != m) {
            throw std::runtime_error(
                    "BinaryEffectiveSourceArchive: m mismatch in " + path.string()
            );
        }

        PatchSide patch_file;
        switch (patch_raw) {
            case 0u:
                patch_file = PatchSide::Left;
                break;
            case 1u:
                patch_file = PatchSide::Right;
                break;
            default:
                throw std::runtime_error(
                        "BinaryEffectiveSourceArchive: bad patch id in " + path.string()
                );
        }

        if (patch_file != patch) {
            throw std::runtime_error(
                    "BinaryEffectiveSourceArchive: patch mismatch in " + path.string()
            );
        }

        const std::size_t nx = static_cast<std::size_t>(nx_u64);
        const std::size_t ny = static_cast<std::size_t>(ny_u64);

        if (nx == 0 || ny == 0) {
            throw std::runtime_error(
                    "BinaryEffectiveSourceArchive: empty axis in " + path.string()
            );
        }

        const std::vector<double> Xd  = read_pod_array<double>(in, nx, "X");
        const std::vector<double> Yd  = read_pod_array<double>(in, ny, "Y");
        const std::vector<double> buf = read_pod_array<double>(in, 2 * nx * ny, "complex data");

        LoadedComponent out;
        out.m = m_file;
        out.patch = patch_file;
        out.component = c_file;
        out.meta.rp = Real(rp_d);
        out.meta.B  = Real(B_d);
        out.meta.sym = symmetry_from_fields_(flags, y_storage_raw);

        out.X.values.resize(nx);
        out.Y.values.resize(ny);
        for (std::size_t i = 0; i < nx; ++i) out.X.values[i] = Real(Xd[i]);
        for (std::size_t j = 0; j < ny; ++j) out.Y.values[j] = Real(Yd[j]);

        out.data.resize(nx * ny);
        for (std::size_t k = 0; k < nx * ny; ++k) {
            out.data[k] = Complex(Real(buf[2 * k]), Real(buf[2 * k + 1]));
        }

        return out;
    }

    std::string BinaryEffectiveSourceArchive::component_tag_(Component c) {
        switch (c) {
            case Component::ll:  return "ll";
            case Component::ln:  return "ln";
            case Component::lm:  return "lm";
            case Component::lmb: return "lmb";
            case Component::nn:  return "nn";
            case Component::nm:  return "nm";
            case Component::nmb: return "nmb";
            case Component::mm:  return "mm";
            case Component::mmb: return "mmb";
            default:
                throw std::runtime_error("BinaryEffectiveSourceArchive: unknown component");
        }
    }

    void BinaryEffectiveSourceArchive::ensure_same_metadata_(const SourceMetadata& a,
                                                             const SourceMetadata& b,
                                                             const std::string& where) {
        if (!same_metadata(a, b)) {
            throw std::runtime_error(where + ": metadata mismatch across files");
        }
    }

    void BinaryEffectiveSourceArchive::ensure_same_axes_(const Axis& a,
                                                         const Axis& b,
                                                         const std::string& where) {
        if (!axis_equal(a, b)) {
            throw std::runtime_error(where);
        }
    }


    SourceSymmetry BinaryEffectiveSourceArchive::symmetry_from_fields_(std::uint32_t flags,
                                                                       std::uint32_t y_storage_raw) {
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
                throw std::runtime_error("BinaryEffectiveSourceArchive: unsupported y_storage id");
        }

        return sym;
    }

} // namespace ghz::source