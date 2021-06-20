#include <iostream>
#include <string>
#include <math.h>
#include <algorithm>
#include "vector.hpp"
#include "scene.hpp"
#include "sphere.hpp"
#include "ray.hpp"
#include "geometry.hpp"
#include "mesh.hpp"

#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

int main() {
	
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();


	int W = 512;
	int H = 512;
	Vector Q(0,0,55);
	Vector S(-10,20,40);
	double alpha = M_PI/3;
	double intensity = 2*pow(10,10);
	double gamma = 2.2;

	std::vector<Geometry*> G;
	//Sphere s1(Vector(0,0,0),10,Vector(1,1,1),false,true,1.5);
	Sphere s2(Vector(0,1000,0), 940, Vector(1,0,0),false,false,1.5);
	Sphere s3(Vector(0,-1000,0), 990, Vector(0,0,1),false,false,1.5);
	Sphere s4(Vector(0,0,-1000), 940, Vector(0,1,0),false,false,1.5);
	Sphere s7(Vector(0,0,1000), 940, Vector(1,0,1),false,false,1.5);
	Sphere s5(Vector(-1000,0,0),940,Vector(0,1,1),false,false,1.5);
	Sphere s6(Vector(1000,0,0),940,Vector(1,1,0),false,false,1.5);
    //Sphere s9(Vector(-20,0,0),10,Vector(1,1,1),true,false,1.5);
	//Sphere s8(Vector(20,0,0),10,Vector(1,1,1),false,true ,1.5);
	//Sphere s10(Vector(20,0,0),9.5,Vector(1,1,1),false,true,1);

	//scene1.pushback(s1);
	G.push_back(new Sphere(Vector(0,1000,0), 940, Vector(1,0,0),false,false,1.5));
	G.push_back(new Sphere(Vector(0,-1000,0), 990, Vector(0,0,1),false,false,1.5));
	G.push_back(new Sphere(Vector(0,0,-1000), 940, Vector(0,1,0),false,false,1.5));
	G.push_back(new Sphere(Vector(0,0,1000), 940, Vector(1,0,1),false,false,1.5));
	G.push_back(new Sphere(Vector(-1000,0,0),940,Vector(0,1,1),false,false,1.5));
	G.push_back(new Sphere(Vector(1000,0,0),940,Vector(1,1,0),false,false,1.5));
	//scene1.pushback(s8);
	//scene1.pushback(s9);
	//scene1.pushback(s10);
	const char cat[14] = "model/cat.obj";
    TriangleMesh* mesh = new TriangleMesh(cat);
	for (int i = 0; i < (mesh->vertices).size(); i++){ 
		(mesh->vertices)[i] = scalar_product((mesh->vertices)[i], 0.6) + Vector(0,-10,0);
	}
	int start = 0;
    int end = (mesh->indices).size();
	mesh->bbox = mesh->compute_box(start, end);
	mesh->compute_BHV(mesh->root, start, end);
	G.push_back(mesh);
	Scene scene1(G, Q, intensity, S);

	int ray_depth = 5;
	double n= 1;
	double K = 1; 
	std::vector<unsigned char> image(W*H * 3, 0);
	#pragma omp parallel for
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			//Ray R = launch_ray(Q, alpha, W, H, i, j);
			//intersection res = intersection_Ray_Scene(R,scene1);
			//Vector L = getColor(res, R, scene1, ray_depth, n);
			Vector L(0,0,0);
			for (int k = 0; k < K; k++){
				Ray R = launch_ray(Q, alpha, W, H, i, j);
				intersection res = scene1.intersect(R);
				L += scene1.getColor(res, R, ray_depth, n);
              
            }
			
			L[0] = L[0] /K;	 
			L[1] = L[1] /K;
			L[2] = L[2] /K;

			double L0_gamma = pow(L[0], 1./gamma) ;
			double L1_gamma = pow(L[1], 1./gamma) ;
			double L2_gamma = pow(L[2], 1./gamma) ;
	
			image[(i*W + j) * 3 + 0] = std::max(std::min(L0_gamma,255.),0.);
			image[(i*W + j) * 3 + 1] = std::max(std::min(L1_gamma,255.),0.);
			image[(i*W + j) * 3 + 2] = std::max(std::min(L2_gamma,255.),0.);
			} 
		
	}
	stbi_write_png("image.png", W, H, 3, &image[0], 0);

	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count() << " milliseconds" << std::endl;

	return 0;
}