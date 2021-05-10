#pragma once
#include <math.h>
#include "vector.hpp"
#include "sphere.hpp"
#include <random>

class Ray {
    public:
        explicit Ray(Vector O, Vector u){
            Origin = O;
            direction = u;
        }
        Vector get_Ray_Origin(){
            return Origin;
        }
        
        Vector get_Ray_direction(){
            return direction;
        }
    private:
        Vector Origin;
        Vector direction;
};

static std::default_random_engine engine(10);
static std::uniform_real_distribution<double> uniform(0,1);
void boxMuller(double stdev, double &x, double &y, double &z){
    double r1 = uniform(engine);
    double r2 = uniform(engine);
    x = sqrt(-2 * log(r1)) * cos(2 * M_PI * r2) * stdev;
    y = sqrt(-2 * log(r1)) * sin(2 * M_PI * r2) * stdev; 
    z = sqrt(r2);
}
Ray launch_ray(const Vector& Q, double alpha, int W, int H, int i, int j){
    int x = j;
    int y = H-i-1;
    double e = 0;
    double f = 0;
    double z = 0; 
    boxMuller(0.5,e,f,z);

    Vector V(x+e+0.5-(W/2),y+f+0.5-(H/2),-W/(2*tan(alpha/2)));  
    Vector V_normalized = vector_normalize(V);

    return Ray(Q,V_normalized);
}

struct intersection {
    bool found = false;
    Vector P;
    Vector N;
    int index;
};

double intersect_ray_plane(Ray& R, Vector A, Vector N){
    
    double uN = dot(R.get_Ray_direction(),N); 
    double t = dot(A-R.get_Ray_Origin(),N)/uN;
    return t;

    
}