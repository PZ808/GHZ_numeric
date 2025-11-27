cmake_minimum_required(VERSION 3.25)
project(GHZ_numeric)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)  # use standard flags, not compiler-specific
# Homebrew default path on Apple Silicon

find_package(OpenMP COMPONENTS CXX)
set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include")
set(OpenMP_CXX_LIBRARIES "-L/opt/homebrew/opt/libomp/lib -lomp")


add_compile_definitions(TEUK_PRECISION=4)
find_package(Eigen3 REQUIRED)
find_package(FFTW3 REQUIRED)


include_directories(
        /opt/homebrew/include
        /opt/homebrew/include/eigen3
)
link_directories(
        /opt/homebrew/lib
)


add_executable(GHZ_numeric main.cpp include/GhzTypes.hpp include/EigenTeukNumTraits.hpp  include/KerrParams.hpp src/KerrParams.cpp
        include/Metric.hpp
        include/Metric.hpp src/KerrMetric.cpp include/KerrMetric.hpp
        include/KerrMetricBL.hpp src/KerrMetricBL.cpp include/MathMacros.hpp include/Coords.hpp src/Coords.cpp #include/Vectors.hpp src/Vectors.cpp
        src/KerrMetricIngoing.cpp include/KerrMetricIngoing.hpp src/KerrMetricOutgoing.cpp
        include/KerrMetricOutgoing.hpp include/Tetrads.hpp
        src/KinnersleyTetrad.cpp include/KinnersleyTetrad.hpp include/SpinCoeffsNP.hpp src/SpinCoeffNP.cpp
        include/WeylScalars.hpp src/WeylScalars.cpp
        include/GHPScalars.hpp src/GHPScalars.cpp include/VectorsGHZ.hpp
        include/KerrMetric.hpp
       # include/SpectralSolverS2.hpp src/SpectralSolverS2.cpp
        include/SpectralSolverS2boost.hpp src/SpectralSolverS2boost.cpp
         include/EllipticIntegrals.hpp  include/KerrOrbit.hpp src/KerrOrbit.cpp
        include/SpectralGHPField.hpp src/SpectralDiffer.cpp src/HeldScalars.cpp include/KerrCharts.hpp)

target_compile_options(GHZ_numeric PRIVATE ${OpenMP_CXX_FLAGS})
target_link_libraries(GHZ_numeric PRIVATE ${OpenMP_CXX_LIBRARIES})

target_link_libraries(GHZ_numeric PRIVATE fftw3l fftw3 m)
