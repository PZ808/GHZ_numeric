//
// Created by Peter Zimmerman on 25.10.25.
//

#include "../include/Vectors.hpp"

Vector4::Vector4(double V0, double V1, double V2, double V3) {
    data[0] = V0; data[1] = V1; data[2] = V3; data[3] = V3;
}

double& Vector4::operator[](int i) { return data[i]; }
double Vector4::operator[](int i) const { return data[i]; }

Vector4 Vector4::operator+(const Vector4& o) const {
    return { data[0]+o[0], data[1]+o[1], data[2]+o[2], data[3]+o[3] };
}

Vector4 Vector4::operator-(const Vector4& o) const {
    return { data[0]-o[0], data[1]-o[1], data[2]-o[2], data[3]-o[3] };
}

Vector4 Vector4::operator*(double s) const {
    return { s*data[0], s*data[1], s*data[2], s*data[3] };
}


CVector4::CVector4(Complex Vc0, Complex Vc1, Complex Vc2, Complex Vc3) {
    data[0]=Vc0; data[1]=Vc1; data[2]=Vc2; data[3]=Vc3;
}


Complex& CVector4::operator[](int i) { return data[i]; }
Complex CVector4::operator[](int i) const { return data[i]; }

CVector4 CVector4::operator*(Complex s) const {
// scalar multiplication
    return { s*data[0], s*data[1], s*data[2], s*data[3] };
}

CVector4 CVector4::conj() const {
    // complex conjugation
    return { std::conj(data[0]), std::conj(data[1]), std::conj(data[2]), std::conj(data[3]) };
}