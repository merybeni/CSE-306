#pragma once
#include <string>
#include <iostream>
#include <stdio.h>
#include <algorithm>
#include <vector>
#include "vector.hpp"
#include "ray.hpp"

class Geometry {
    public:
        Vector color;
        bool mirror;
        bool transparent; 
        double indice_refrac;

        virtual intersection intersect(Ray& R) = 0;
};

