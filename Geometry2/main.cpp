#include "svg.cpp"
#include "sample.cpp"




std::vector<Vector> generate_points(int n) {
	std::vector<Vector> points(n);
	for (int i = 0; i < n; i++) {
		points[n] = Vector((double) rand() / RAND_MAX, (double) rand() / RAND_MAX);
	}
	return points;
}

Polygon generate_polygon(int n) {
	Polygon Points;
	for (int i = 0; i < n; i++) {
		Points.add( Vector((double) rand() / RAND_MAX, (double) rand() / RAND_MAX));
	}
	return Points;
}

std::vector<double> generate_weights(int n) {
	std::vector<double> weights(n);
	for (int i = 0; i < n; i++) {
		weights[i] = (double) rand() / RAND_MAX;
	}
	return weights;
}

std::vector<double> normal_lambda(std::vector<Vector> points) {
    Vector center(0.5, 0.5);
    std::vector<double> lambdas;
    for (int i = 0; i < int(points.size()); i++) {
        double lambda = exp(-pow(norm(points[i] - center), 2) / 0.02);
        lambdas.push_back(lambda / double(points.size()));
    }
    return lambdas;
}

int main() {
    Polygon subjectPolygon;
    subjectPolygon.add(Vector(0.7,0.2,0));
    subjectPolygon.add(Vector(0.5,0.2,0));
    subjectPolygon.add(Vector(0.3,0.8,0));
    subjectPolygon.add(Vector(0.5,0.8,0));


    Polygon clipPolygon;
    clipPolygon.add(Vector(0.1,0.6,0));
    clipPolygon.add(Vector(0.9,0.6,0));
    clipPolygon.add(Vector(0.9,0.4,0));
    clipPolygon.add(Vector(0.1,0.4,0));

    std::vector<Polygon> polygons;
    Polygon clippedPolygon = SutherlandHodgman(subjectPolygon,clipPolygon);
    //polygons.push_back(subjectPolygon);
    //polygons.push_back(clipPolygon);
    polygons.push_back(clippedPolygon);
    save_svg(polygons,"poly.svg");


    Polygon bounds;
	bounds.vertices.push_back(Vector(0,0,0));
	bounds.vertices.push_back(Vector(0,1,0));
	bounds.vertices.push_back(Vector(1,1,0));
	bounds.vertices.push_back(Vector(1,0,0));

  
    Polygon Points = generate_polygon(1000);
    save_svg(Voronoi(bounds, Points), "voronoi.svg");
	std::vector<double> weights = generate_weights(1000);
	save_svg(PVoronoi(bounds, Points, weights), "powervoronoi.svg");
   

	std::vector<Vector> points2 = generate_points(90);
    Polygon Points2 = generate_polygon(90);
	std::vector<double> lambdas = normal_lambda(points2);
	std::vector<double> weights2 = OptimalTransport(Points2, lambdas);
	save_svg(PVoronoi(bounds, Points2, weights2), "optimaltransport.svg");

    return 0;
}






