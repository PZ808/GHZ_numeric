//
// Created by Peter Zimmerman on 27.10.25.
//
#include "../include/TeukTypes.hpp"
#include "../include/VectorsGHZ.hpp"
#include "../include/Tetrads.hpp"


void TetradTransformations::conformal_transform(Tetrad &tetrad, const Real Om) {

    auto common = Om;
    auto l_tilde =  ghz::Vector4{tetrad.l[0], tetrad.l[1], tetrad.l[2], tetrad.l[3]} * (1.0/(Om*Om));
    auto n_tilde =  ghz::Vector4{tetrad.n[0], tetrad.n[1], tetrad.n[2], tetrad.n[3]} * Om*Om;

    tetrad.l = l_tilde;
    tetrad.n = n_tilde;
}