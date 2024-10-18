// loop-subdivision
//
// Author   : Omkar Vengurlekar (Arizona State University)
// Email    : ovengurl@asu.edu

#include <fstream>
#include <vector>
#include <cmath>
#include <unordered_map>

#include "OBJFileReader.h"
#include "Solid.h"
#include "iterators.h"
#include "SolidDelegate.h"
 
using namespace MeshLib;


int main(int argc, char *argv[])
{
	// Read in the obj file
	Solid mesh;
	OBJFileReader of;
	std::ifstream in(argv[1]);
	of.readToSolid(&mesh, in);

	/******************* Put you subdivision processing here *********************/
	Solid newMesh;
	int vid = 0;
	int fid = 0;
	SolidDelegate delegate;

	std::unordered_map<Vertex*, Vertex*> vertexMap;
	for ( SolidVertexIterator viter(&mesh); !viter.end(); ++viter){
		Vertex *v = *viter;
		Vertex *vNew = delegate.createVertex(&newMesh, ++vid);

		if (v->boundary()){
			std::vector <Point > plist;
			HalfEdge *hedge = v->halfedge();

			if (hedge->edge()->boundary()){
				do {
					hedge = hedge->clw_rotate_about_target();
				} while (hedge->clw_rotate_about_target());
				if (hedge == hedge->edge()->halfedge(0)){
					hedge = hedge->edge()->halfedge(1);
				}
				else if (hedge == hedge->edge()->halfedge(1)){
					hedge = hedge->edge()->halfedge(0);
				}
				hedge = hedge->clw_rotate_about_source();
				vNew->point() = v->point()*3/4 + (v->halfedge()->source()->point() + hedge->target()->point())/8;
			}
			else{
				do { 
					hedge = hedge->ccw_rotate_about_target();
					plist.push_back(hedge->source()->point()); 
				} while (!hedge->edge()->boundary());
				vNew->point() = v->point() * 3.0f / 4.0f + plist.back() / 8.0f;
				do { 
					hedge = hedge->clw_rotate_about_target();
					plist.push_back(hedge->source()->point()); 
				} while (hedge->clw_rotate_about_target());
				if (hedge == hedge->edge()->halfedge(0)){
					hedge = hedge->edge()->halfedge(1);
				}
				else if (hedge == hedge->edge()->halfedge(1)){
					hedge = hedge->edge()->halfedge(0);
				}
				hedge = hedge->clw_rotate_about_source();
				vNew->point() = vNew->point() + hedge->target()->point() / 8.0f;
			}
		}
		else{
			std::vector <Point > plist;
			HalfEdge *hedge = v->halfedge();

			int n = 0;
			do {
				plist.push_back(hedge->source()->point());
				hedge = hedge->clw_rotate_about_target();
				n++;
			} while (hedge != v->halfedge());

			float alpha;
			if (n > 3){
				float center = (0.375f + (0.25f * cos((2.0f * 3.1415926f) / (float)n)));
				alpha = (0.625f - (center * center)) / (float)n;
			}
			else {
				alpha = 3.0f / 16.0f;
			}

			Point temp = { 0.0f, 0.0f, 0.0f };
			for (int i = 0; i < n; i++){
				temp += plist.back();
				plist.pop_back();
			}
			vNew->point() = v->point() *(1 - n*alpha) + temp.operator*(alpha);
		}
		vertexMap[v] = vNew;
	}

	std::unordered_map<Edge*, Vertex*> edgeMap;
	SolidEdgeIterator eiter(&mesh);
	for (; !eiter.end(); ++eiter){
		Edge *e = *eiter;
		Vertex *ev1 = mesh.idVertex(e->vertex(0));
		Vertex *ev2 = mesh.idVertex(e->vertex(1));
		Vertex *vNew = delegate.createVertex(&newMesh, ++vid);

		if (e->boundary()){
			vNew->point() = (ev1->point() + ev2->point()) / 2.0f;
		}
		else{
			vNew->point() = (ev1->point() + ev2->point()) * 3.0f / 8.0f;
			Vertex *vt1 = e->halfedge(0)->he_next()->target();
			Vertex *vt2 = e->halfedge(1)->he_next()->target();
			vNew->point() = vNew->point() + (vt1->point() + vt2->point()) / 8.0f;
		}
		edgeMap[e] = vNew;
	} 

	SolidFaceIterator faceiter(&mesh);
	for(; !faceiter.end(); ++faceiter){
		Face *f = *faceiter;
		HalfEdge * fhe[3];
		fhe[0] = f->halfedge();
		fhe[1] = fhe[0]->he_next();
		fhe[2] = fhe[1]->he_next();

		int vIndex[3];
		for (int i = 0; i < 3; i++){
			vIndex[i] = edgeMap[fhe[i]->edge()]->id();
		}
		delegate.createFace(&newMesh, vIndex, vid);

		for (int i = 0; i < 3; i++){
			vIndex[0] = vertexMap[fhe[i]->source()]->id();
			vIndex[1] = edgeMap[fhe[i]->edge()]->id();
			vIndex[2] = edgeMap[fhe[(i+2)%3]->edge()]->id();
			delegate.createFace(&newMesh, vIndex, ++fid);
		}
	}

	newMesh.labelBoundaryEdges(); 

	// Write out the resultant obj file
	int vObjID = 1;
	std::map<int, int> vidToObjID;

	std::ofstream os(argv[2]);

	SolidVertexIterator iter(&newMesh);

	for(; !iter.end(); ++iter)
	{
		Vertex *v = *iter;
		Point p = v->point();
		os << "v " << p[0] << " " << p[1] << " " << p[2] << std::endl;
		vidToObjID[v->id()] = vObjID++;
	}
	os << "# " << (unsigned int)newMesh.numVertices() << " vertices" << std::endl;

	float u = 0.0, v = 0.0;
	for(iter.reset(); !iter.end(); ++iter)
	{
		Vertex *vv = *iter;
		std::string key( "uv" );
		std::string s = Trait::getTraitValue (vv->string(), key );
		if( s.length() > 0 )
		{
			sscanf( s.c_str (), "%f %f", &u, &v );
		}
		os << "vt " << u << " " << v << std::endl;
	}
	os << "# " << (unsigned int)newMesh.numVertices() << " texture coordinates" << std::endl;

	SolidFaceIterator fiter(&newMesh);
	for(; !fiter.end(); ++fiter)
	{
		Face *f = *fiter;
		FaceVertexIterator viter(f);
		os << "f " ;
		for(; !viter.end(); ++viter)
		{
			Vertex *v = *viter;
			os << vidToObjID[v->id()] << "/" << vidToObjID[v->id()] << " ";
		}
		os << std::endl;
	}
	os.close();

}