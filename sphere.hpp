#pragma once
#include <math.h>
#include "vector.hpp"
#include "ray.hpp"
#include "geometry.hpp"

class Sphere : public Geometry {
    public:
        explicit Sphere(Vector C, double R,Vector albedo, bool M, bool T, double n){
            Center = C;
            Radius = R;
            color = albedo;
            mirror = M;
            transparent = T;
            indice_refrac = n; 
        }

        explicit Sphere(){}

        Vector get_Sphere_Center(){
            return Center;
        }
        double get_Sphere_Radius(){
            return Radius;
        }
      
        intersection intersect(Ray& R) {
            intersection res;
            Vector u = R.get_Ray_direction();
            Vector O = R.get_Ray_Origin();
            Vector C = Center; 
            double radius = Radius;
            Vector OC = O-C;

            double discriminant = pow(dot(u,OC),2) - pow(norm(OC),2)+ pow(radius,2);

            if (discriminant == 0) {
                double t = dot(u,C-O);
                if (t >= 0) {
                    res.P = O + scalar_product(u,t);
                    res.N = scalar_product(res.P-C,1/norm(res.P-C));
                    res.found = true;
                }
            }
            if (discriminant > 0) {
                double t1 = dot(u,C-O) - sqrt(discriminant);
                double t2 = dot(u,C-O) + sqrt(discriminant);
                if (t1 >= 0){
                    res.P = O + scalar_product(u,t1);
                    res.N = scalar_product(res.P-C,1/norm(res.P-C));
                    res.found = true;
                } else {
                    if (t1<=0 && t2 >0) {
                        res.P = O + scalar_product(u,t2);
                    res.N = scalar_product(res.P-C,1/norm(res.P-C));
                    res.found = true;
                    }
                }
            }
            return res;   
        }


    private:
        Vector Center;
        double Radius;
};


Vector random_cos(const Vector& N){
    double r1 = uniform(engine);
    double r2 = uniform(engine);
    double x = cos(2*M_PI*r1) * sqrt(1-r2);
    double y = sin(2*M_PI*r1) * sqrt(1-r2);
    double z = sqrt(r2);
    double min = INFINITY;
    double index;
    for (int i = 0; i < 3; i++) {
        if (abs(N[i]) < min) { 
            min = abs(N[i]);
            index = i;
        }
    }
    Vector T1;
    if (index == 0) {
        T1 = Vector(0, N[2], -N[1]);
    } else {
        if (index == 1) {
            T1 = Vector(N[2], 0, -N[0]);
        } else {
            T1 = Vector(N[1], -N[0], 0);
        }
    }
    T1 = vector_normalize(T1);
    Vector T2 = cross(N,T1);
    return scalar_product(T1,x) + scalar_product(T2,y) + scalar_product(N,z);
}