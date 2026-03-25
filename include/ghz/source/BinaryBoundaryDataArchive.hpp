//
// Created by Peter Zimmerman on 17.03.26.
//
#ifndef GHZ_SOURCE_BINARYBOUNDARYDATAARCHIVE_HPP
#define GHZ_SOURCE_BINARYBOUNDARYDATAARCHIVE_HPP

#pragma once

#include "ghz/source/EffectiveSource.hpp"

#include <array>
#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

namespace ghz::source {

    enum class BoundaryComponent : std::uint32_t {
        delta = 0,
        deltaPrime = 1,
        count
    };

    inline constexpr std::size_t boundary_component_count =
            static_cast<std::size_t>(BoundaryComponent::count);

    // loaded dataset for one mode and one patch (left or right)
    struct BoundaryMetadata {
        teuk::Real r_boundary{0};
        teuk::Real rp{0}; // radius of particle (circular orbit)
        teuk::Real B{0};  // B converts from X to rp
        SourceSymmetry sym{}; // symmetry of the source mode, see SourceSymmetry enum
    };

    struct BoundaryModeData {
        int m{0};                         // azimuthal mode number
        PatchSide patch{PatchSide::Left}; // patch side (left or right)
        Axis Y;                           // Y-axis values (angular coordinate)


        std::array<std::vector<teuk::Complex>, boundary_component_count> components{};

        [[nodiscard]] std::vector<teuk::Complex>& component(BoundaryComponent c) {
            return components[static_cast<std::size_t>(c)];
        }

        [[nodiscard]] const std::vector<teuk::Complex>& component(BoundaryComponent c) const {
            return components[static_cast<std::size_t>(c)];
        }
    };

    class BinaryBoundaryDataArchive {
    public:
        using Real = teuk::Real;
        using Complex = teuk::Complex;

        BinaryBoundaryDataArchive(std::filesystem::path root_dir,
                                  std::string prefix,
                                  int probe_m,
                                  PatchSide probe_patch);

        [[nodiscard]] const BoundaryMetadata& metadata() const noexcept { return metadata_; }

        [[nodiscard]] BoundaryModeData load_mode_patch(int m, PatchSide patch) const;

        [[nodiscard]] BoundaryModeData load_nonnegative_mode_patch_raw_(int m,
                                                                        PatchSide patch) const;

    private:
        struct LoadedBoundaryComponent {
            int m{0};
            PatchSide patch{PatchSide::Left};
            BoundaryComponent component{BoundaryComponent::delta};

            BoundaryMetadata meta;
            Axis Y;
            std::vector<Complex> data;
        };

        std::filesystem::path root_dir_;
        std::string prefix_;
        BoundaryMetadata metadata_{};

        [[nodiscard]] std::filesystem::path component_path_(BoundaryComponent c,
                                                            int m,
                                                            PatchSide patch) const;

        [[nodiscard]] LoadedBoundaryComponent read_component_file_(BoundaryComponent c,
                                                                   int m,
                                                                   PatchSide patch) const;

        static std::string component_tag_(BoundaryComponent c);

        static void ensure_same_metadata_(const BoundaryMetadata& a,
                                          const BoundaryMetadata& b,
                                          const std::string& where);

        static void ensure_same_axes_(const Axis& a,
                                      const Axis& b,
                                      const std::string& where);

        static SourceSymmetry symmetry_from_fields_(std::uint32_t flags,
                                                    std::uint32_t y_storage_raw);
    };

} // namespace ghz::source

#endif // GHZ_SOURCE_BINARYBOUNDARYDATAARCHIVE_HPP