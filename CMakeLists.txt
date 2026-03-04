cmake_minimum_required(VERSION 3.20)
project(GHZ_numeric LANGUAGES CXX)

# --------------------------
# Compiler settings
# --------------------------
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Force Homebrew LLVM
set(CMAKE_C_COMPILER   "/opt/homebrew/opt/llvm@20/bin/clang")
set(CMAKE_CXX_COMPILER "/opt/homebrew/opt/llvm@20/bin/clang++")

add_compile_definitions(TEUK_PRECISION=2)

# --------------------------
# Homebrew include and link paths
# --------------------------
include_directories(
        /opt/homebrew/include
        /opt/homebrew/include/eigen3
        /opt/homebrew/opt/libomp/include  # OpenMP headers
)
link_directories(
        /opt/homebrew/lib
        /opt/homebrew/opt/libomp/lib      # OpenMP library
)

# --------------------------
# OpenMP (Homebrew LLVM)
# --------------------------
set(OpenMP_CXX_FLAGS "-fopenmp -I/opt/homebrew/opt/libomp/include")
set(OpenMP_C_FLAGS   "-fopenmp -I/opt/homebrew/opt/libomp/include")
set(OpenMP_CXX_LIB_NAMES "omp")
set(OpenMP_C_LIB_NAMES "omp")
set(OpenMP_omp_LIBRARY "/opt/homebrew/opt/libomp/lib/libomp.dylib")

find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    message(STATUS "Found OpenMP: ${OpenMP_CXX_FLAGS}")
endif()

# --------------------------
# External dependencies
# --------------------------
find_package(Eigen3 REQUIRED)
find_package(FFTW3 REQUIRED)

# --------------------------
# Executable
# --------------------------
add_executable(GHZ_numeric
        main.cpp

        include/GhzTypes.hpp
        include/EigenTeukNumTraits.hpp
        include/KerrParams.hpp            src/KerrParams.cpp
        include/Metric.hpp                src/KerrMetric.cpp
        include/KerrMetric.hpp
        include/KerrMetricBL.hpp          src/KerrMetricBL.cpp
        include/MathMacros.hpp
        include/Coords.hpp                src/Coords.cpp
        include/KerrMetricIngoing.hpp     src/KerrMetricIngoing.cpp
        include/KerrMetricOutgoing.hpp    src/KerrMetricOutgoing.cpp

        include/Tetrads.hpp               src/KinnersleyTetrad.cpp
        include/KinnersleyTetrad.hpp

        include/SpinCoeffsNP.hpp          src/SpinCoeffNP.cpp
        include/WeylScalars.hpp           src/WeylScalars.cpp
        include/GHPScalars.hpp            src/GHPScalars.cpp

        include/VectorsGHZ.hpp

        include/SpectralSolverS2boost.hpp src/SpectralSolverS2boost.cpp

        include/EllipticIntegrals.hpp
        include/KerrOrbit.hpp             src/KerrOrbit.cpp

        sand/SpectralGHPField.hpp
        include/SpectralDiffer.hpp        src/SpectralDiffer.cpp
        include/HeldScalars.hpp           src/HeldScalars.cpp
        include/KerrCharts.hpp

        include/SpectralGHPFieldVectorized.hpp
        include/Splines.hpp src/KinnersleySpectralHeldOperators.cpp
        include/ODE.hpp
        include/Corrector.hpp)

# --------------------------
# Target include directories
# --------------------------
target_include_directories(GHZ_numeric PRIVATE include)

# --------------------------
# Link libraries
# --------------------------
target_link_libraries(GHZ_numeric PRIVATE
        OpenMP::OpenMP_CXX
        fftw3 fftw3l m
        )
