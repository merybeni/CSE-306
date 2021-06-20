#pragma once
#include <iostream>
#include <math.h>
#include <vector>
#include <random>
#include "vector.hpp"
#include "sphere.hpp"
#include "ray.hpp"

class Scene {
    public:
        explicit Scene(std::vector<Geometry*> init, Vector center_camera, double intensity, Vector source){
            scene = init;
            Q = center_camera;
            I = intensity;
            S = source;
        }
       
        Vector get_Q() {
            return Q;
        }
        std::vector<Geometry*> get_Scene(){
            return scene;
        }
        Vector get_ligh_source(){
            return S;
        }
        double get_intensity(){
            return I;
        }
        intersection intersect(Ray& R) {
            intersection result;
            double smallest_distance = INFINITY;
            for (int i = 0; i < scene.size(); i++) {
                intersection res = scene[i]->intersect(R);
                double D = distance(res.P, R.get_Ray_Origin());
                if (res.found == true && D < smallest_distance) {
                    res.index = i;
                    result = res;
                    smallest_distance = D;
                }
            }
            return result;
        }
        Vector compute_shadow(intersection& result, Ray R){
            Vector N = result.N;
            Vector P = result.P;
            Geometry* Sphere = scene[result.index];
            Vector albe = Sphere->color;
           
            Vector D = S-P;
            double d = norm(D);
            Vector w_i = vector_normalize(D);
            double V_P;
            Vector P_ep = P + scalar_product(N,0.00001);
            Ray ray_from_P = Ray(P_ep,w_i);
            intersection visib = this->intersect(ray_from_P);
           

            if (visib.found == true){
                Vector P_1 = visib.P;
                double t = norm(visib.P-P);
            
                if (t > d) {
                   
                    V_P = 1;
                }
                else {
                    V_P = 0;
                }
            }
            if (visib.found == false) {
                V_P = 1;
            }

            

            double e_1 = (4* pow(M_PI,2) * (pow(d,2)));
            double j = 0;
            double e_2 = I * V_P * std::max(dot(N,w_i),j);

            double e_3 = e_2/e_1;
            Vector L = scalar_product(albe,e_3);
            return L;
            
        }
        Vector getColor(intersection result, Ray& ray, int ray_depth, double n) {
            Vector w_i = ray.get_Ray_direction();
            Vector w_r = w_i;
        
            if (ray_depth < 0){
                return Vector(0,0,0);
            }
            Geometry* sphere_id = scene[result.index];
            Vector P = result.P;
            Vector N = result.N;
            Vector P_ep = P + scalar_product(N,0.00001);
            double n_sphere = sphere_id->indice_refrac;
            Vector albedo =  sphere_id->color;
                
            
            if (result.found == true){
                

                if (sphere_id->mirror == true){
                    
                    Vector w_r = w_i - scalar_product(N,2*dot(w_i,N));
                    Vector P_ep = P + scalar_product(N,0.00001);
                    Ray reflected_ray = Ray(P_ep,w_r);
                    result = this->intersect(reflected_ray);
                    return getColor(result,reflected_ray, ray_depth-1,n);
                }
                
                else if  (sphere_id->transparent == true){
                    
                    double n1;
                    double n2;

                    

                    if (dot(N,w_r)< 0){
                        n1 = n;
                        n2 = n_sphere; 
                        n = n2;

                    }

                    else {
                        n1 = n;
                        n2 = 1;
                        n = n2; 
                        N = scalar_product(N,-1);
                    }

                    double k0 = pow(n1-n2,2) / pow(n1+n2,2);
                    double R = k0 + (1-k0)* pow(1 - abs(dot(N,w_r)),5);
                    double u = (double) rand()/RAND_MAX;    
                    
                    if (u<R){
                        Vector w_r = w_i - scalar_product(N,2*dot(w_i,N));
                        Vector P_ep = P + scalar_product(N,0.00001);
                        Ray reflected_ray = Ray(P_ep,w_r);
                        result = this->intersect(reflected_ray);
                        return getColor(result,reflected_ray, ray_depth-1,n);
                        }

                    else {

                        double sin = sqrt(1-pow(dot(w_r,N),2));  
                        if (n1 > n2 && sin > (n2/n1)){
                            Vector w_r = w_i - scalar_product(N,2*dot(w_i,N));
                            Vector P_ep = P + scalar_product(N,0.00001);
                            Ray reflected_ray = Ray(P_ep,w_r);
                            result = this->intersect(reflected_ray);
                            return getColor(result,reflected_ray, ray_depth-1,n);  
                        }
                        else {
                            Vector w_rt = scalar_product(w_r -  scalar_product(N,dot(w_r,N)), n1/n2);
                            Vector w_rn =  scalar_product( N,-1 * sqrt(1- pow(n1/n2,2) * (1- pow(dot(w_r,N),2))));
                            w_r = w_rt + w_rn;

                            Vector P_ep = P - scalar_product(N,0.00001);
                            Ray refracted_ray = Ray(P_ep,w_r);
                            result = this->intersect(refracted_ray);
                            return getColor(result,refracted_ray, ray_depth-1,n);
                        }
                    }
                }

                else {
                    
                    Vector L0  =  this->compute_shadow(result, ray); 
                    
                
                    Vector randomvect = random_cos(N);
                    Ray randomRay = Ray(P_ep, randomvect); 
                    result = this->intersect(randomRay);
                
                    L0 += vector_multiplication(albedo, getColor(result, randomRay, ray_depth-1, n ));
                    return L0;
                }
            }

            else {
                return Vector(0,0,0);
            }
        }

   


    private:
        std::vector<Geometry*> scene;
        Vector Q;
        Vector S;
        double I;
};