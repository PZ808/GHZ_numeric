//
// Created by Peter Zimmerman on 31.10.25.
//

#ifndef GHZ_NUMERIC_CORRECTOR_HPP
#define GHZ_NUMERIC_CORRECTOR_HPP

class HeldScalar;
struct HeldCoefficients;
class Trajectory;
class Source;

class Corrector:
private:
    GHPScalar aH_, bH_, cH_, eH_, fH_;
public:
    Corrector(HeldCoefficients hc, Trajectory traj, Source src);

#endif //GHZ_NUMERIC_CORRECTOR_HPP
