#pragma once
#include <string>
#include <iostream>
#include <stdio.h>
#include <algorithm>
#include <vector>
#include <list>
#include "vector.hpp"
#include "ray.hpp"
#include "geometry.hpp"

class BoundingBox{
    public: 
        Vector Bmin; 
        Vector Bmax;

    public: 
        BoundingBox(){};

        
        bool intersect(Ray& ray, double distance){

            Vector T0;
            Vector T1;

            Vector O= ray.get_Ray_Origin();
            Vector u= ray.get_Ray_direction();
            for (int i = 0; i < 3; i++){
                T0[i] = (Bmin[i]-O[i])/u[i];
                T1[i] = (Bmax[i]-O[i])/u[i];
          
            }
            double t0_x,t0_y,t0_z,t1_x,t1_y,t1_z;
		    
            if(T0[0]<T1[0]) {
                
                t0_x = T0[0];
                t1_x=T1[0];} 

            else{
                
                t0_x = T1[0];
                t1_x=T0[0];}
            
            if(T0[1]<T1[1]) {
                
                t0_y = T0[1];
                t1_y=T1[1];} 
            
            else{
                
                t0_y = T1[1];
                t1_y=T0[1];}
            
            if(T0[2]<T1[2]) {
                
                t0_z = T0[2];
                t1_z=T1[2];} 
            
            else{t0_z = T1[2];
                
                t1_z=T0[2];}
		    

		
		
		    double l = std::max(std::max(t0_x,t0_y),t0_z);
		    if (std::min(std::min(t1_x,t1_y),t1_z)>l){
                distance = l;
                return true;
            }
            else{
                return false;}
            
            
        }
    
};

class Node {
    public:
   
        Node* left;
        Node* right;
        BoundingBox bbox;
        int starting_triangle;
        int ending_triangle;
     

    Node() {
        this->left = NULL;
        this->right= NULL;
    }
    
};

class TriangleIndices {
public:

    TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1, int group = -1, bool added = false) : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk), group(group) {
    };
    int vtxi, vtxj, vtxk; 
    int uvi, uvj, uvk;  
    int ni, nj, nk;  
    int group;       
};
 
 
class TriangleMesh : public Geometry {

public:
    std::vector<TriangleIndices> indices;
    std::vector<Vector> vertices;
    std::vector<Vector> normals;
    std::vector<Vector> uvs;
    std::vector<Vector> vertexcolors;
    BoundingBox bbox;
    Node* root;

    ~TriangleMesh() {}

    TriangleMesh(const char* obj) {
        color = Vector(1,1,1);
        mirror = false; 
        transparent = false; 
        indice_refrac = 1; 
        this->readOBJ(obj);
        this->root = new Node();

    }

    BoundingBox compute_box(int& starting_triangle, int& ending_triangle){
        BoundingBox Box1;
        Box1.Bmin = Vector(INFINITY, INFINITY,INFINITY);
        Box1.Bmax = Vector( -INFINITY, -INFINITY, -INFINITY);
        for (int i = starting_triangle; i < ending_triangle; i++){
            TriangleIndices I = indices[i];
            Vector A  = vertices[I.vtxi];
            Vector B  = vertices[I.vtxj];
            Vector C  = vertices[I.vtxk];

            for (int j = 0; j < 3; j++){
                if (std::min(std::min(A[j],B[j]),C[j]) < Box1.Bmin[j]){
                    Box1.Bmin[j] = std::min(std::min(A[j],B[j]),C[j]);
                }
                if (std::max(std::max(A[j],B[j]),C[j]) > Box1.Bmax[j]){
                    Box1.Bmax[j] = std::max(std::max(A[j],B[j]),C[j]);
                }
            }
        }
        return Box1;
    }

    Vector compute_barycenter(int& i) {
        TriangleIndices T = indices[i];
        return scalar_product(vertices[T.vtxi] + vertices[T.vtxj] + vertices[T.vtxk], pow(3,-1));
    }

    int get_longest(Vector& diag) {
        int index = 0;
        double max = diag[0];
        for (int i = 1; i < 3; i++){ 
            if (diag[i] > max){ 
                max = diag[i];
                index = i;
            }
        }   
        return index;
    }

    void compute_BHV(Node* node, int& start, int& end) {
        node->bbox = this->compute_box(start, end);
        node->starting_triangle = start;
        node->ending_triangle = end;
        Vector diag = node->bbox.Bmax - node->bbox.Bmin;
        Vector middle_diag = node->bbox.Bmin + scalar_product(diag, 0.5);
        int longest_axis = get_longest(diag);
        int pivot_index = start;
        for (int i = start ; i < end ; i++) {
            Vector barycenter = this->compute_barycenter(i);
            if (barycenter[longest_axis] < middle_diag[longest_axis]) { 
                std::swap(indices[i], indices[pivot_index]); 
                pivot_index++;
            }
        }
        if (pivot_index <= start || pivot_index >= end-1 || end-start < 5) {
            return;
        }
        node->left = new Node();
        node->right = new Node();
        compute_BHV(node->left, start, pivot_index); 
        compute_BHV(node->right, pivot_index, end);
    }

    virtual intersection intersect(Ray& ray) { 
        intersection res;
        double nothing;
        if (!root->bbox.intersect(ray, nothing)) {
            return res;
        }
        std::list<Node*> nodes_to_visit; 
        nodes_to_visit.push_front(root);
        double best_inter_distance = INFINITY;
        intersection I;
        Node* curNode;
        while (!nodes_to_visit.empty()) {
            curNode = nodes_to_visit.back();
            nodes_to_visit.pop_back();
            if (!curNode->left) {
                double inter_distance;
                if (curNode->left->bbox.intersect(ray, inter_distance)) {
                    if(inter_distance < best_inter_distance) {
                        nodes_to_visit.push_back(curNode->left);
                    }
                }
                if (curNode->right->bbox.intersect(ray, inter_distance)) {
                    if(inter_distance < best_inter_distance) {
                        nodes_to_visit.push_back(curNode->right);
                    }
                } 
            } else {
                double masafa = INFINITY; 
                for (int i = 0; i < indices.size(); i++){
                    TriangleIndices I = indices[i];
                    Vector A = vertices[I.vtxi];
                    Vector B = vertices[I.vtxj];
                    Vector C = vertices[I.vtxk];
                    Vector u = ray.get_Ray_direction();
                    Vector O = ray.get_Ray_Origin();
                    Vector e1 = B - A; 
                    Vector e2 = C - A; 
                    Vector N = cross(e1,e2);
                    double uN = dot(u,N);
                    double beta = dot(e2,cross((A-O),u))/ uN;
                    double gamma = - dot(e1,cross((A-O),u))/ uN;
                    double alpha = 1 - beta - gamma; 
                    if (uN != 0){
                        if( beta>0 && gamma >0 && alpha >0){
                            double t = dot(A-O,N)/dot(u,N);
                            if (t < masafa && t > 0) {
                                masafa = t;
                                res.found = true; 
                                res.P = O + scalar_product(u,t); 
                                res.N = vector_normalize(N);
                            }
                        }
                    }
                }
                if (masafa < best_inter_distance){
                    best_inter_distance = masafa;
                    I = res;
                }
            }
        }
        return I;
    }
    
    
    
    void readOBJ(const char* obj) {
 
        char matfile[255];
        char grp[255];
 
        FILE* f;
        f = fopen(obj, "r");
        int curGroup = -1;
        while (!feof(f)) {
            char line[255];
            if (!fgets(line, 255, f)) break;
 
            std::string linetrim(line);
            linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
            strcpy(line, linetrim.c_str());
 
            if (line[0] == 'u' && line[1] == 's') {
                sscanf(line, "usemtl %[^\n]\n", grp);
                curGroup++;
            }
 
            if (line[0] == 'v' && line[1] == ' ') {
                Vector vec;
 
                Vector col;
                if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec[0], &vec[1], &vec[2], &col[0], &col[1], &col[2]) == 6) {
                    col[0] = std::min(1., std::max(0., col[0]));
                    col[1] = std::min(1., std::max(0., col[1]));
                    col[2] = std::min(1., std::max(0., col[2]));
 
                    vertices.push_back(vec);
                    vertexcolors.push_back(col);
 
                } else {
                    sscanf(line, "v %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
                    vertices.push_back(vec);
                }
            }
            if (line[0] == 'v' && line[1] == 'n') {
                Vector vec;
                sscanf(line, "vn %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
                normals.push_back(vec);
            }
            if (line[0] == 'v' && line[1] == 't') {
                Vector vec;
                sscanf(line, "vt %lf %lf\n", &vec[0], &vec[1]);
                uvs.push_back(vec);
            }
            if (line[0] == 'f') {
                TriangleIndices t;
                int i0, i1, i2, i3;
                int j0, j1, j2, j3;
                int k0, k1, k2, k3;
                int nn;
                t.group = curGroup;
 
                char* consumedline = line + 1;
                int offset;
 
                nn = sscanf(consumedline, "%u/%u/%u %u/%u/%u %u/%u/%u%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
                if (nn == 9) {
                    if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                    if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                    if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                    if (j0 < 0) t.uvi = uvs.size() + j0; else   t.uvi = j0 - 1;
                    if (j1 < 0) t.uvj = uvs.size() + j1; else   t.uvj = j1 - 1;
                    if (j2 < 0) t.uvk = uvs.size() + j2; else   t.uvk = j2 - 1;
                    if (k0 < 0) t.ni = normals.size() + k0; else    t.ni = k0 - 1;
                    if (k1 < 0) t.nj = normals.size() + k1; else    t.nj = k1 - 1;
                    if (k2 < 0) t.nk = normals.size() + k2; else    t.nk = k2 - 1;
                    indices.push_back(t);
                } else {
                    nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
                    if (nn == 6) {
                        if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                        if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                        if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                        if (j0 < 0) t.uvi = uvs.size() + j0; else   t.uvi = j0 - 1;
                        if (j1 < 0) t.uvj = uvs.size() + j1; else   t.uvj = j1 - 1;
                        if (j2 < 0) t.uvk = uvs.size() + j2; else   t.uvk = j2 - 1;
                        indices.push_back(t);
                    } else {
                        nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2, &offset);
                        if (nn == 3) {
                            if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                            if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                            if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                            indices.push_back(t);
                        } else {
                            nn = sscanf(consumedline, "%u//%u %u//%u %u//%u%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
                            if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                            if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                            if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                            if (k0 < 0) t.ni = normals.size() + k0; else    t.ni = k0 - 1;
                            if (k1 < 0) t.nj = normals.size() + k1; else    t.nj = k1 - 1;
                            if (k2 < 0) t.nk = normals.size() + k2; else    t.nk = k2 - 1;
                            indices.push_back(t);
                        }
                    }
                }
 
                consumedline = consumedline + offset;
 
                while (true) {
                    if (consumedline[0] == '\n') break;
                    if (consumedline[0] == '\0') break;
                    nn = sscanf(consumedline, "%u/%u/%u%n", &i3, &j3, &k3, &offset);
                    TriangleIndices t2;
                    t2.group = curGroup;
                    if (nn == 3) {
                        if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                        if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                        if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                        if (j0 < 0) t2.uvi = uvs.size() + j0; else  t2.uvi = j0 - 1;
                        if (j2 < 0) t2.uvj = uvs.size() + j2; else  t2.uvj = j2 - 1;
                        if (j3 < 0) t2.uvk = uvs.size() + j3; else  t2.uvk = j3 - 1;
                        if (k0 < 0) t2.ni = normals.size() + k0; else   t2.ni = k0 - 1;
                        if (k2 < 0) t2.nj = normals.size() + k2; else   t2.nj = k2 - 1;
                        if (k3 < 0) t2.nk = normals.size() + k3; else   t2.nk = k3 - 1;
                        indices.push_back(t2);
                        consumedline = consumedline + offset;
                        i2 = i3;
                        j2 = j3;
                        k2 = k3;
                    } else {
                        nn = sscanf(consumedline, "%u/%u%n", &i3, &j3, &offset);
                        if (nn == 2) {
                            if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                            if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                            if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                            if (j0 < 0) t2.uvi = uvs.size() + j0; else  t2.uvi = j0 - 1;
                            if (j2 < 0) t2.uvj = uvs.size() + j2; else  t2.uvj = j2 - 1;
                            if (j3 < 0) t2.uvk = uvs.size() + j3; else  t2.uvk = j3 - 1;
                            consumedline = consumedline + offset;
                            i2 = i3;
                            j2 = j3;
                            indices.push_back(t2);
                        } else {
                            nn = sscanf(consumedline, "%u//%u%n", &i3, &k3, &offset);
                            if (nn == 2) {
                                if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                                if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                                if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                                if (k0 < 0) t2.ni = normals.size() + k0; else   t2.ni = k0 - 1;
                                if (k2 < 0) t2.nj = normals.size() + k2; else   t2.nj = k2 - 1;
                                if (k3 < 0) t2.nk = normals.size() + k3; else   t2.nk = k3 - 1;                             
                                consumedline = consumedline + offset;
                                i2 = i3;
                                k2 = k3;
                                indices.push_back(t2);
                            } else {
                                nn = sscanf(consumedline, "%u%n", &i3, &offset);
                                if (nn == 1) {
                                    if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                                    if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                                    if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                                    consumedline = consumedline + offset;
                                    i2 = i3;
                                    indices.push_back(t2);
                                } else {
                                    consumedline = consumedline + 1;
                                }
                            }
                        }
                    }
                }
 
            }
 
        }
        fclose(f);
 
    }
 
    
    
};