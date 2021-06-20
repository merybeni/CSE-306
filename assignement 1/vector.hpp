#pragma once
#include <string>
#include <math.h>
#include <algorithm>

class Vector {
    public:
        explicit Vector(double x=0, double y=0, double z=0){
            coords[0] = x;
            coords[1] = y;
            coords[2] = z;
        }
        Vector& operator += (const Vector& b){
            coords[0] += b[0];
            coords[1] += b[1];
            coords[2] += b[2];  
            return *this;    
        }
        const double& operator[](int i) const {return coords[i];}
        double& operator[](int i) {return coords[i];}
    private:
        double coords[3];
};   
    
Vector operator+(const Vector& a,const Vector& b){
    return Vector(a[0]+b[0],a[1]+b[1],a[2]+b[2]);
}

Vector operator-(const Vector& a,const Vector& b){
    return Vector(a[0]-b[0],a[1]-b[1],a[2]-b[2]);
}

double dot(const Vector& a,const Vector& b){
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

Vector cross(const Vector& a,const Vector& b){
    return Vector(a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]);
}

double norm(const Vector& a){
    return sqrt(pow(a[0],2)+pow(a[1],2)+pow(a[2],2));
}

Vector vector_normalize(const Vector& a){
    double n = norm(a);
    if (!n) {
        return a;
    }
    return Vector(a[0]/n,a[1]/n,a[2]/n);
}

Vector scalar_product(const Vector& a, double s){
    return Vector(s*a[0],s*a[1],s*a[2]);
}

Vector get_vector_two_points(const Vector& a,const Vector& b){
    return Vector(b[0]-a[0],b[1]-a[1],b[2]-a[2]);
}

double distance (const Vector& a,const Vector& b){
    return sqrt(pow(a[0]-b[0],2)+pow(a[1]-b[1],2)+pow(a[2]-b[2],2));
}

Vector vector_multiplication(const Vector& a,const Vector& b){
    return Vector(a[0]*b[0], a[1]*b[1], a[2]*b[2]);
}

