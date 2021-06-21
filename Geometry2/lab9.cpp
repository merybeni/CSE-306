#include <set>
#include <vector>
#include "vector.hpp"
#include  "mesh.hpp"
#include <cmath>



std::vector<Vector> Tutte(TriangleMesh T, int n_iter){
    int n_size = T.vertices.size();
    int T_indices = T.indices.size();

    std::vector<Vector> delta_M(n_size); // identify ordered boundary vertices 
    // boundary length 

    double sum = 0;
    for (int i=0; i< n_size; i++){
        sum += norm(delta_M[i+1]-delta_M[i]);
    }

    double cs =0; 
    std::vector<Vector> VO(T_indices);
    std::vector<double> theta(n_size); 

    for (int i; i< n_size; i++){
        theta[i]= 2 * M_PI * (cs/sum); 
        VO[i] = Vector(cos(theta[i]), sin(theta[i]), 0);
        cs+= norm(delta_M[i+1]-delta_M[i]);
        
    }

    for (int iter=0; iter< n_iter; iter++){
        for (int i= n_size; i< T_indices; i++){
            Vector adj_indice(T.indices[i].vtxi, T.indices[i].vtxj, T.indices[i].vtxk);
            std::vector<Vector> adj_i(T_indices);
            for (int j = 0; j<T_indices; j++){
                
                for (int k=0; k<3; k++){
                    adj_i[adj_indice[k]].insert(T.vertices[(k+1)%3]);
                    adj_i[adj_indice[k]].insert(T.vertices[(k+2)%3]);    
                }

                for (int h=0; h<adj_i.size(); h++){
                    VO[i] += adj_i[h];
                }

                VO[i] = (1/adj_i.size()) *  VO[i]

                

             }

        }

    }


    return VO[n_iter-1];

}


//T.vtxi T.vtxj, T.vtk 
