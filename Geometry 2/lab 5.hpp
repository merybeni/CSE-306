#include <math.h>
#include <random>
#include "vector.hpp"
#include <iostream>
#include "time.h"
#include <iomanip>

#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

static std::default_random_engine engine(10);
static std::uniform_real_distribution<double> uniform(0,1);
void  Global_Illumination( double &x, double &y, double &z){
    double r1 = uniform(engine);
    double r2 = uniform(engine);
    x = sqrt(r2*(1-r2)) * cos(2 * M_PI * r1) ;
    y = sqrt(r2*(1-r2)) * sin(2 * M_PI * r1) ; 
    z = 1-2*r2;
}

void color_matching(unsigned char* input_image, unsigned char* output_image, int niter, int input_W, int input_H, int input_C, int model_C) {
    int n = input_W * input_H;
    for (int i = 0; i < niter; i++) {
        double x;
        double y;
        double z; 
        Global_Illumination(x, y, z);
        Vector V(x,y,z);
        
        //We store the dot product and pixel index as a pair of values
        std::vector<std::pair<double, int> > projsource(n);
        std::vector<std::pair<double, int> > projtarget(n);

        for (int j=0; j<n ; j++){
            Vector point_source(input_image[j * input_C],input_image[j * input_C + 1],input_image[j * input_C + 2]);
            Vector point_target(output_image[j * model_C],output_image[j * model_C +1],output_image[j * model_C + 2]);

            projsource[j] =  std::make_pair(dot(point_source,V),j);
            projtarget[j] =  std::make_pair(dot(point_target,V),j);
        }

        // Sort according to the dot product 
        std::sort(projsource.begin(),projsource.end());
        sort(projtarget.begin(),projsource.end());

        // Advect initial point cloud.
        for (int i=0; i<n ; i++) {
            Vector L(input_image[projsource[i].second * input_C], input_image[projsource[i].second * input_C + 1],input_image[projsource[i].second * input_C + 2]);
            Vector d((projtarget[i].first - projsource[i].first)*V[0], (projtarget[i].first - projsource[i].first)*V[1] ,(projtarget[i].first - projsource[i].first)*V[2]);
            add_update(L,d);
            input_image[projsource[i].second*input_C] = L[0];
            input_image[projsource[i].second*input_C + 1] = L[1];
            input_image[projsource[i].second*input_C + 2] = L[2];
        }
    }
    stbi_write_png("output.png", input_W, input_H, input_C, &input_image[0], 0);
}

int main() {
    int input_W, input_H, input_C;
    unsigned char *input_image = stbi_load("input.png", &input_W, &input_H, &input_C, 0);

    int model_W, model_H, model_C;
    unsigned char *output_image = stbi_load("model.png", &model_W, &model_H, &model_C, 0);

    int niter = 200;

    color_matching(input_image, output_image, niter, input_W, input_H, input_C, model_C);
    return 0;
}