#include <stdio.h>
#include "lbfgs.h"
#include "svg.cpp"
#include <cmath>


double integral ( Polygon Poly , Vector P_i){
    // Poly.vertices = {(X_i, Y_i)}_i 
    // Poly.vertices[i] = (X_i, Y_i)   

    double integral = 0;
    for (int i = 0; i < Poly.vertices.size(); i ++) {
        double Xk = Poly.vertices[i][0];
        double Yk = Poly.vertices[i][1];
        double Xk_1 = Poly.vertices[i-1][0];
        double Yk_1 = Poly.vertices[i-1][1];
        integral += ( Xk_1* Yk - Xk * Yk_1) * (Xk_1*Xk_1 + Xk_1*Xk + pow(Xk,2) + Yk_1*Yk_1 + Yk_1*Yk + pow(Yk,2)
                - 4 * (P_i[0]*(Xk_1 + Xk) + P_i[1]*(Yk_1 + Yk)) + 6 * pow(norm(P_i),2));
    }
    return integral /12; 
}

double function_g(Polygon subjectpolygon, Polygon Points, std::vector<double> weights, std::vector<double> lambdas){
    


    std::vector<Polygon> polygons= PVoronoi(subjectpolygon, Points, weights);
    double sum = 0;
    for (int i=0; i<Points.vertices.size();i++){
        
         sum += integral(polygons[i], Points.vertices[i]) - weights[i] * (polygons[i].area()) + lambdas[i] * weights[i];

    }
    return sum;

}

double gradient(std::vector<Polygon> polygons,  std::vector<double> lambdas, int i){
    double grad = -  polygons[i].area() + lambdas[i];
    return grad;

}

class objective_function
{
protected:
    lbfgsfloatval_t *m_x;
    Polygon Points;
    std::vector<double> lambdas;

public:
    objective_function() : m_x(NULL)
    {
    }

    virtual ~objective_function()
    {
        if (m_x != NULL) {
            lbfgs_free(m_x);
            m_x = NULL;
        }
    }

    std::vector<double> run(Polygon polygon, std::vector<double> lambdas)
    {   
        
        lbfgsfloatval_t fx;
        double N = lambdas.size();
        lbfgsfloatval_t *m_x = lbfgs_malloc(N);
        this->lambdas = lambdas;
        this->Points = polygon;


        /* Initialize the variables. */
        for (int i = 0;i < N;i += 2) {
            m_x[i] = 0.001;
        
        }

        /*
            Start the L-BFGS optimization; this will invoke the callback functions
            evaluate() and progress() when necessary.
         */
        int ret = lbfgs(N, m_x, &fx, _evaluate, NULL, this, NULL);

        /* Report the result. */
        printf("L-BFGS optimization terminated with status code = %d\n", ret);
        printf("  fx = %f, x[0] = %f, x[1] = %f\n", fx, m_x[0], m_x[1]);
        
        std::vector<double> result(N);
        for (int i=0; i<N; i++){
            result[i] = m_x[i];
        }

        return result;
        
    }

protected:
    static lbfgsfloatval_t _evaluate(
        void *instance,
        const lbfgsfloatval_t *x,
        lbfgsfloatval_t *g,
        const int n,
        const lbfgsfloatval_t step
        )
    {
        return reinterpret_cast<objective_function*>(instance)->evaluate(x, g, n, step);
    }

    lbfgsfloatval_t evaluate(
        const lbfgsfloatval_t *x,
        lbfgsfloatval_t *g,
        const int n,
        const lbfgsfloatval_t step
        )
    {
        lbfgsfloatval_t fx = 0.0;

        std::vector<double> weights(n);
        
        for (int i = 0; i < n; i ++) {
            weights[i] = x[i];
        }
        
        Polygon subjectpolygon;
        subjectpolygon.vertices.push_back(Vector(0,0,0));
        subjectpolygon.vertices.push_back(Vector(0,1,0));
        subjectpolygon.vertices.push_back(Vector(1,1,0));
        subjectpolygon.vertices.push_back(Vector(1,0,0));

        std::vector<Polygon> polygons = PVoronoi(subjectpolygon, Points, weights);
        
        for (int i = 0; i < n ; i ++) {
            g[i] = - gradient(polygons, lambdas,i);
        }

        fx = - function_g(subjectpolygon, Points, weights, lambdas);

        return fx;
    }
};

std::vector<double> OptimalTransport(Polygon Points, std::vector<double> lambdas) {
    objective_function obj;
    return obj.run(Points, lambdas);
}
