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
    
Vector operator+(const Vector& a,const Vector &b){
    return Vector(a[0]+b[0],a[1]+b[1],a[2]+b[2]);
}

Vector operator-(const Vector& a,const Vector &b){
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

void add_update( Vector& a, Vector& b){
    a[0] = a[0] + b[0];
    a[1] = a[1] + b[1];
    a[2] = a[2] + b[2];
    
}

Vector intersect(std::vector<Vector> edge, Vector A, Vector B ){
    Vector V1 = edge[0];
    Vector V2 = edge[1];
    Vector N(V2[1]-V1[1],V1[0]-V2[0],0);

    double t = dot(V1-A,N)/dot(B-A,N);
    Vector P;
    if (t<0 || t>1){
        return P;
    }
    P = A+scalar_product(B-A,t);
    return P;
}


bool inside(std::vector<Vector> edge, Vector P ){
    Vector V1 = edge[0];
    Vector V2 = edge[1];
    Vector N;
    N[0] = V2[1]-V1[1];
    N[1] = V1[0]-V2[0];

    if (dot(P-V1,N)<= 0){
        return true;
    }
    return false;
}


Vector intersectVoronoi(Vector A, Vector B, Vector Pi, Vector Pj) {
    Vector M = scalar_product(Pi+Pj,0.5);
    Vector P;
    double t = dot(M-A, Pi-Pj) / dot(B-A, Pi-Pj);
    P = A + scalar_product(B-A,t );
    return P;
}

bool insideVoronoi(Vector X, Vector Pi, Vector Pj) {
    Vector M = scalar_product(Pi+Pj, 0.5);
    if (dot(X-M, Pj-Pi) < 0) {
        return true;
    }
    return false;
}


Vector intersectpowerVoronoi(Vector V1, Vector V2, Vector Pi, Vector Pj, double wi, double wj){
  
    Vector M = scalar_product(Pi + Pj,0.5) + scalar_product(Pj-Pi,(wi-wj)/(2*(pow(norm(Pi-Pj),2))));

    double t = dot(M-V1,Pi-Pj)/dot(V2-V1,Pi-Pj);
    Vector P = V1 +scalar_product(V2-V1,t);
    return P;
}

bool insidepowerVoronoi(Vector X, Vector Pi, Vector Pj, double wi, double wj ){
  
    Vector M = scalar_product(Pi + Pj,0.5) + scalar_product(Pj-Pi,(wi-wj)/(2*(pow(norm(Pi-Pj),2))));

    if (dot(X-M,Pj-Pi)< 0){
        return true;
    }
    return false;
}


