#ifndef GHZ_SOURCE_BINARY_EFFECTIVE_SOURCE_ARCHIVE_HPP
#define GHZ_SOURCE_BINARY_EFFECTIVE_SOURCE_ARCHIVE_HPP
#pragma once

#include "ghz/source/EffectiveSource.hpp"

#include <filesystem>
#include <string>

namespace ghz::source {

    class BinaryEffectiveSourceArchive final : public EffectiveSourceArchive {
    public:
        BinaryEffectiveSourceArchive(std::filesystem::path root_dir,
                                     std::string prefix = "Teff",
                                     int probe_m = 0,
                                     PatchSide probe_patch = PatchSide::Left);

        SourceMetadata metadata() const override { return metadata_; }

    protected:
        ModeData load_nonnegative_mode_patch_raw_(int m, PatchSide patch) const override;

    private:
        struct LoadedComponent {
            int m = 0;
            PatchSide patch = PatchSide::Left;
            Component component = Component::ll;
            SourceMetadata meta{};
            Axis X;
            Axis Y;
            std::vector<Complex> data;
        };

        std::filesystem::path root_dir_;
        std::string prefix_;
        SourceMetadata metadata_{};

        std::filesystem::path component_path_(Component c, int m, PatchSide patch) const;
        LoadedComponent read_component_file_(Component c, int m, PatchSide patch) const;

        static std::string component_tag_(Component c);

        static void ensure_same_metadata_(const SourceMetadata& a,
                                          const SourceMetadata& b,
                                          const std::string& where);

        static void ensure_same_axes_(const Axis& a,
                                      const Axis& b,
                                      const std::string& where);

        static SourceSymmetry symmetry_from_fields_(std::uint32_t flags,
                                                    std::uint32_t y_storage_raw);
    };

} // namespace ghz::source

#endif
