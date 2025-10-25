//
// Created by Peter Zimmerman on 25.10.25.
//

#ifndef GHZ_NUMERIC_VECTORS_HPP
#define GHZ_NUMERIC_VECTORS_HPP
#pragma once
#include <complex>
#include <cmath>


using Complex = std::complex<double>;


struct Vector4 {
    double data[4];
    Vector4(double V0=0, double V1=0, double V2=0, double V3=0);


    double& operator[](int i);
    double operator[](int i) const;


    Vector4 operator+(const Vector4& o) const;
    Vector4 operator-(const Vector4& o) const;
    Vector4 operator*(double s) const;
};


struct CVector4 {
    Complex data[4];
    CVector4(Complex V0=0, Complex V1=0, Complex V2=0, Complex V3=0);


    Complex& operator[](int i);
    Complex operator[](int i) const;


    CVector4 operator*(Complex s) const;
    CVector4 conj() const;
};

#endif //GHZ_NUMERIC_VECTORS_HPP
