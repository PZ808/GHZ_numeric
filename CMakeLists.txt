cmake_minimum_required(VERSION 3.25)
project(GHZ_numeric)

set(CMAKE_CXX_STANDARD 17)
find_package(Eigen3 REQUIRED)
find_package(FFTW3 REQUIRED)

include_directories(
        /opt/homebrew/include
        /opt/homebrew/include/eigen3
)

link_directories(
        /opt/homebrew/lib
)

add_executable(GHZ_numeric main.cpp include/KerrParams.hpp src/KerrParams.cpp include/Metric.hpp include/Metric.hpp src/KerrMetric.cpp include/KerrMetric.hpp
        include/KerrMetricBL.hpp src/KerrMetricBL.cpp include/MathMacros.hpp include/Coords.hpp src/Coords.cpp #include/Vectors.hpp src/Vectors.cpp
        src/KerrMetricIngoing.cpp include/KerrMetricIngoing.hpp src/KerrMetricOutgoing.cpp include/KerrMetricOutgoing.hpp include/Tetrads.hpp
        src/KinnersleyTetrad.cpp include/KinnersleyTetrad.hpp include/SpinCoeffsNP.hpp src/SpinCoeffNP.cpp include/WeylScalars.hpp src/WeylScalars.cpp
        include/GHPScalars.hpp src/GHPScalars.cpp include/TeukTypes.hpp include/VectorsGHZ.hpp include/SpectralSolverS2.hpp src/SpectralSolverS2.cpp)


target_link_libraries(GHZ_numeric PRIVATE fftw3 m)

