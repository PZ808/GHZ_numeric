cmake_minimum_required(VERSION 3.25)
project(GHZ_numeric)

set(CMAKE_CXX_STANDARD 17)

add_executable(GHZ_numeric main.cpp include/KerrParams.hpp src/KerrParams.cpp include/Metric.hpp include/Metric.hpp src/KerrMetric.cpp include/KerrMetric.hpp
        include/KerrMetricBL.hpp src/KerrMetricBL.cpp include/MathMacros.hpp include/Coords.hpp src/Coords.cpp include/Vectors.hpp src/Vectors.cpp
        src/KerrMetricIngoing.cpp include/KerrMetricIngoing.hpp src/KerrMetricOutgoing.cpp include/KerrMetricOutgoing.hpp include/Tetrads.hpp
        src/KinnersleyTetrad.cpp include/KinnersleyTetrad.hpp include/SpinCoeffsNP.hpp src/SpinCoeffNP.cpp)
        #Coords.cpp Coords.hpp KerrMetric.cpp KerrMetric.hpp Metric.hpp HeldScalars.cpp HeldScalars.hpp GHPScalars.cpp GHPScalars.hpp MathMacros.hpp)
