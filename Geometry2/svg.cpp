#pragma once
#include <iostream>
#include <vector>
#include "vector.hpp"

// if the Polygon class name conflicts with a class in wingdi.h on Windows, use a namespace or change the name
class Polygon {  
    public:
        std::vector<Vector> vertices;
		Polygon(){};
        Polygon(std::vector<Vector> V){
            vertices = V;
        }

        void add(Vector v){
            vertices.push_back(v);
        }

		double area(){
			double air = 0;
			for (int i=0; i<vertices.size(); i++){
				if (i==0){
					air += (vertices[vertices.size()][0] + vertices[0][1]) * (vertices[vertices.size()][1] * vertices[0][1]);
				}
				air += (vertices[i-1][0] + vertices[i][0]) * (vertices[i-1][1] - vertices[i][1]);
			}
			return abs(air / 2);
}


		
		
};


Polygon SutherlandHodgman(Polygon subjectpolygon, Polygon clippolygon){
    Polygon OutPolygon;
    std::vector<std::vector<Vector> > clipedges;
    for (int i=0; i<clippolygon.vertices.size(); i++) {
        Vector V1 = clippolygon.vertices[i];
        Vector V2 = clippolygon.vertices[(i > 0) ? (i - 1) : (clippolygon.vertices.size() - 1)];
        std::vector<Vector> clipedge;
        clipedge.push_back(V1);
        clipedge.push_back(V2);
        clipedges.push_back(clipedge);
    }

    for (int i=0; i<clippolygon.vertices.size();i++){
		OutPolygon = Polygon();
        for (int j=0; j<subjectpolygon.vertices.size(); j++){
            Vector curVertex = subjectpolygon.vertices[j];
            Vector prevVertex = subjectpolygon.vertices[(j > 0)?(j-1):(subjectpolygon.vertices.size()-1)];
            
            Vector intersection = intersect(clipedges[i], prevVertex, curVertex);
            if (inside(clipedges[i], curVertex )== true){
                if (inside(clipedges[i], prevVertex )== false){
                    OutPolygon.add(intersection);   
                }
                OutPolygon.add(curVertex); 
            }
            else if(inside(clipedges[i], prevVertex )== true){
                OutPolygon.add(intersection);
            }

        }
        subjectpolygon = OutPolygon;
    }
    return OutPolygon;
}




std::vector<Polygon> Voronoi(Polygon subjectpolygon, Polygon Points){
	std::vector<Polygon> Voronoicells;
	Polygon cell;
	for (int i=0; i<Points.vertices.size(); i++){
		cell = subjectpolygon;
		Vector Point_i = Points.vertices[i];
		
		for (int j=0; j<Points.vertices.size(); j++){
			if (j!=i){
				Polygon OutPolygon = Polygon();
				Vector Point_j = Points.vertices[j];
				for (int z=0; z<cell.vertices.size(); z++){
					Vector curVertex = cell.vertices[z];
            		Vector prevVertex = cell.vertices[(z > 0)?(z-1):(cell.vertices.size()-1)];
            
            		Vector intersection = intersectVoronoi(prevVertex, curVertex, Point_i, Point_j);
					if (insideVoronoi(curVertex,Point_i,Point_j)== true){
						if (insideVoronoi(prevVertex, Point_i, Point_j)==false){
							OutPolygon.add(intersection);
						}
						OutPolygon.add(curVertex);
					}

					else if (insideVoronoi(prevVertex, Point_i, Point_j)==true){
						OutPolygon.add(intersection);
					}
				}
				
				cell = OutPolygon; 
			}
		}
		Voronoicells.push_back(cell);
		
	}
	return Voronoicells;
}


std::vector<Polygon> PVoronoi(Polygon subjectpolygon, Polygon Points, std::vector<double> weights){
	std::vector<Polygon> Voronoicells;
	Polygon cell;
	for (int i=0; i<Points.vertices.size(); i++){
		cell = subjectpolygon;
		Vector Point_i = Points.vertices[i];
		double weight_i = weights[i];
		
		for (int j=0; j<Points.vertices.size(); j++){
			if (j!=i){
				Polygon OutPolygon = Polygon();
				Vector Point_j = Points.vertices[j];
				double weight_j = weights[j];
				for (int z=0; z<cell.vertices.size(); z++){
					Vector curVertex = cell.vertices[z];
            		Vector prevVertex = cell.vertices[(z > 0)?(z-1):(cell.vertices.size()-1)];
            
            		Vector intersection = intersectpowerVoronoi(prevVertex, curVertex, Point_i, Point_j, weight_i, weight_j);
					if (insidepowerVoronoi(curVertex,Point_i,Point_j,weight_i, weight_j)== true){
						if (insidepowerVoronoi(prevVertex, Point_i, Point_j, weight_i, weight_j)==false){
							OutPolygon.add(intersection);
						}
						OutPolygon.add(curVertex);
					}

					else if (insidepowerVoronoi(prevVertex, Point_i, Point_j, weight_i, weight_j)==true){
						OutPolygon.add(intersection);
					}
				}
				
				cell = OutPolygon; 
			}
		}
		Voronoicells.push_back(cell);
		
	}
	return Voronoicells;
}


// saves a static svg file. The polygon vertices are supposed to be in the range [0..1], and a canvas of size 1000x1000 is created
	void save_svg(const std::vector<Polygon> &polygons, std::string filename, std::string fillcol = "none") {
		FILE* f = fopen(filename.c_str(), "w+"); 
		fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\" height = \"1000\">\n");
        for (int i=0; i<polygons.size(); i++) {
		    fprintf(f, "<g>\n");
		    fprintf(f, "<polygon points = \""); 
		    for (int j = 0; j < polygons[i].vertices.size(); j++) {
			    fprintf(f, "%3.3f, %3.3f ", (polygons[i].vertices[j][0] * 1000), (1000 - polygons[i].vertices[j][1] * 1000));
		    }
		    fprintf(f, "\"\nfill = \"%s\" stroke = \"black\"/>\n", fillcol.c_str());
		    fprintf(f, "</g>\n");
        }
		fprintf(f, "</svg>\n");
		fclose(f);
	}


// Adds one frame of an animated svg file. frameid is the frame number (between 0 and nbframes-1).
// polygons is a list of polygons, describing the current frame.
// The polygon vertices are supposed to be in the range [0..1], and a canvas of size 1000x1000 is created
	void save_svg_animated(const std::vector<Polygon> &polygons, std::string filename, int frameid, int nbframes) {
		FILE* f;
		if (frameid == 0) {
			f = fopen(filename.c_str(), "w+");
			fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\" height = \"1000\">\n");
			fprintf(f, "<g>\n");
		} else {
			f = fopen(filename.c_str(), "a+");
		}
		fprintf(f, "<g>\n");
		for (int i = 0; i < polygons.size(); i++) {
			fprintf(f, "<polygon points = \""); 
			for (int j = 0; j < polygons[i].vertices.size(); j++) {
				fprintf(f, "%3.3f, %3.3f ", (polygons[i].vertices[j][0] * 1000), (1000-polygons[i].vertices[j][1] * 1000));
			}
			fprintf(f, "\"\nfill = \"none\" stroke = \"black\"/>\n");
		}
		fprintf(f, "<animate\n");
		fprintf(f, "	id = \"frame%u\"\n", frameid);
		fprintf(f, "	attributeName = \"display\"\n");
		fprintf(f, "	values = \"");
		for (int j = 0; j < nbframes; j++) {
			if (frameid == j) {
				fprintf(f, "inline");
			} else {
				fprintf(f, "none");
			}
			fprintf(f, ";");
		}
		fprintf(f, "none\"\n	keyTimes = \"");
		for (int j = 0; j < nbframes; j++) {
			fprintf(f, "%2.3f", j / (double)(nbframes));
			fprintf(f, ";");
		}
		fprintf(f, "1\"\n	dur = \"5s\"\n");
		fprintf(f, "	begin = \"0s\"\n");
		fprintf(f, "	repeatCount = \"indefinite\"/>\n");
		fprintf(f, "</g>\n");
		if (frameid == nbframes - 1) {
			fprintf(f, "</g>\n");
			fprintf(f, "</svg>\n");
		}
		fclose(f);
	}