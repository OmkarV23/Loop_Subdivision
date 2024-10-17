#include <fstream>
#include <vector>
#include <unordered_map>

#include "OBJFileReader.h"
#include "Solid.h"
#include "SolidDelegate.h"
#include "iterators.h"

using namespace MeshLib;


std::vector<Face*> get_adjacent_faces(Solid *mesh, Edge *edge) {
    std::vector<Face*> adjacent_faces;

    HalfEdge *he1 = edge->halfedge(0);
    if (he1 != nullptr && he1->face() != nullptr) {
        adjacent_faces.push_back(he1->face());
    }
    HalfEdge *he2 = edge->halfedge(1);
    if (he2 != nullptr && he2->face() != nullptr) {
        adjacent_faces.push_back(he2->face());
    }

    return adjacent_faces;
}

std::vector<Vertex*> get_neighbor_vertices(Vertex* v) {
    std::vector<Vertex*> neighbors;

    HalfEdge* he = v->halfedge();
    if (he == nullptr) {
        return neighbors;
    }

    HalfEdge* start_he = he;
    do {
        Vertex* neighbor = he->target();
        neighbors.push_back(neighbor);

        HalfEdge* he_sym = he->he_sym();
        if (he_sym == nullptr) {
            break;
        }

        he = he_sym->he_next();
        if (he == nullptr || he == start_he) {
            break;
        }
    } while (true);

    return neighbors;
}


bool is_boundary_vertex(Vertex* v) {
    HalfEdge* he = v->halfedge();
    if (he == nullptr) {
        return true;
    }

    HalfEdge* start_he = he;
    do {
        if (he->face() == nullptr || he->he_sym() == nullptr) {
            return true;
        }

        HalfEdge* he_sym = he->he_sym();
        if (he_sym == nullptr) {
            return true;
        }
        he = he_sym->he_next();
        if (he == nullptr || he == start_he) {
            break;
        }
    } while (true);

    return false;
}

void get_boundary_neighbors(Vertex* v, Vertex*& v_prev, Vertex*& v_next) {
    v_prev = nullptr;
    v_next = nullptr;

    HalfEdge* he = v->halfedge();
    if (he == nullptr) {
        return; // Isolated vertex
    }
    HalfEdge* boundary_he = nullptr;
    HalfEdge* start_he = he;
    do {
        if (he->he_sym() == nullptr || he->face() == nullptr) {
            boundary_he = he;
            break;
        }

        he = he->he_sym()->he_next();
        if (he == nullptr || he == start_he) {
            break;
        }
    } while (true);

    if (boundary_he == nullptr) {
        return; // No boundary half-edge found
    }
    v_next = boundary_he->target();

    HalfEdge* he_prev = boundary_he->he_prev();
    if (he_prev != nullptr) {
        v_prev = he_prev->source();
    }
}

MeshLib::Edge* getEdgeBetweenVertices(MeshLib::Solid* mesh, MeshLib::Vertex* v1, MeshLib::Vertex* v2) {
    MeshLib::Edge* e = mesh->vertexEdge(v1, v2);
    if (e == nullptr) {
        e = mesh->vertexEdge(v2, v1);
    }
    return e;
}


int main(int argc, char *argv[])
{
	// Read in the obj file
	Solid mesh;
	OBJFileReader of;
	std::ifstream in(argv[1]);
	of.readToSolid(&mesh, in);

	std::cout << "Read in " << mesh.numVertices() << " vertices, " << mesh.numEdges() << " edges, and " << mesh.numFaces() << " faces." << std::endl;

	/******************* Put you subdivision processing here *********************/


    std::unordered_map<MeshLib::Edge*, MeshLib::Vertex*> edgeMidpoints;
    MeshLib::SolidDelegate delegate;

    std::unordered_map<MeshLib::Vertex*, int> vertexToID;
    int nextVertexID = 1;

    for (MeshLib::SolidVertexIterator viter(&mesh); !viter.end(); ++viter) {
        MeshLib::Vertex* v = *viter;
        int id = v->id();
        vertexToID[v] = id;
        nextVertexID = max(nextVertexID, id + 1);
    }

    for (MeshLib::SolidEdgeIterator eiter(&mesh); !eiter.end(); ++eiter) {
        MeshLib::Edge* e = *eiter;
        MeshLib::Vertex* v1 = e->halfedge(0)->source();
        MeshLib::Vertex* v2 = e->halfedge(0)->target();

        MeshLib::Vertex* v3 = nullptr;
        MeshLib::Vertex* v4 = nullptr;

        MeshLib::HalfEdge* he1 = e->halfedge(0);
        if (he1->face() != nullptr) {
            v3 = he1->he_next()->target();
        }

        MeshLib::HalfEdge* he2 = e->halfedge(1);
        if (he2 != nullptr && he2->face() != nullptr) {
            v4 = he2->he_next()->target();
        }

        MeshLib::Point new_pos;
        if (v3 != nullptr && v4 != nullptr) {
            // Internal edge
            new_pos = (v1->point() + v2->point()) * (3.0 / 8.0) +
                      (v3->point() + v4->point()) * (1.0 / 8.0);
        } else {
            // Boundary edge
            new_pos = (v1->point() + v2->point()) * 0.5;
        }

        MeshLib::Vertex* new_vertex = delegate.createVertex(&mesh, nextVertexID++);
        new_vertex->point() = new_pos;
        edgeMidpoints[e] = new_vertex;
        vertexToID[new_vertex] = new_vertex->id();
    }

    std::unordered_map<MeshLib::Vertex*, MeshLib::Point> newVertexPositions;

    for (MeshLib::SolidVertexIterator viter(&mesh); !viter.end(); ++viter) {
        MeshLib::Vertex* v = *viter;
        if (is_boundary_vertex(v)) {
            MeshLib::Vertex* v_prev;
            MeshLib::Vertex* v_next;
            get_boundary_neighbors(v, v_prev, v_next);

            if (v_prev != nullptr && v_next != nullptr) {
                MeshLib::Point new_pos = v->point() * (3.0 / 4.0) +
                                         (v_prev->point() + v_next->point()) * (1.0 / 8.0);
                newVertexPositions[v] = new_pos;
            } else {
                newVertexPositions[v] = v->point();
            }
        } else {
            std::vector<MeshLib::Vertex*> neighbors = get_neighbor_vertices(v);
            int n = neighbors.size();
            double beta;
            if (n == 3) {
                beta = 3.0 / 16.0;
            } else {
                beta = (1.0 / n) * ((5.0 / 8.0) - std::pow((3.0 / 8.0 + 0.25 * std::cos(2.0 * M_PI / n)), 2.0));
            }
            MeshLib::Point sum_neighbors(0.0, 0.0, 0.0);
            for (MeshLib::Vertex* neighbor : neighbors) {
                sum_neighbors += neighbor->point();
            }
            MeshLib::Point new_pos = v->point() * (1.0 - n * beta) + sum_neighbors * beta;
            newVertexPositions[v] = new_pos;
        }
    }


    for (auto& pair : newVertexPositions) {
        MeshLib::Vertex* v = pair.first;
        v->point() = pair.second;
    }

    std::vector<MeshLib::Face*> oldFaces;
    for (MeshLib::SolidFaceIterator fiter(&mesh); !fiter.end(); ++fiter) {
        MeshLib::Face* face = *fiter;
        oldFaces.push_back(face);
    }

    int nextFaceID = 1;
    for (MeshLib::SolidFaceIterator fiter(&mesh); !fiter.end(); ++fiter) {
        MeshLib::Face* f = *fiter;
        int id = f->id();
        nextFaceID = max(nextFaceID, id + 1);
    }

    MeshLib::Solid newMesh;

    for (auto& pair : vertexToID) {
        MeshLib::Vertex* v = pair.first;
        int id = pair.second;
        MeshLib::Vertex* new_v = delegate.createVertex(&newMesh, id);
        new_v->point() = v->point();
        vertexToID[v] = id;
    }

    int faces_created = 0;
    for (MeshLib::Face* face : oldFaces) {
        MeshLib::HalfEdge* he0 = face->halfedge();
        MeshLib::HalfEdge* he1 = he0->he_next();
        MeshLib::HalfEdge* he2 = he1->he_next();

        MeshLib::Vertex* v0 = he0->source();
        MeshLib::Vertex* v1 = he1->source();
        MeshLib::Vertex* v2 = he2->source();

        MeshLib::Edge* e0 = getEdgeBetweenVertices(&mesh, v0, v1);
        MeshLib::Edge* e1 = getEdgeBetweenVertices(&mesh, v1, v2);
        MeshLib::Edge* e2 = getEdgeBetweenVertices(&mesh, v2, v0);

        if (e0 == nullptr || e1 == nullptr || e2 == nullptr) {
            std::cerr << "Error: One or more edges are missing for face ID: " << face->id() << std::endl;
            continue;
        }

        // Check if edge midpoints exist
        if (edgeMidpoints.find(e0) == edgeMidpoints.end() ||
            edgeMidpoints.find(e1) == edgeMidpoints.end() ||
            edgeMidpoints.find(e2) == edgeMidpoints.end()) {
            std::cerr << "Error: Missing edge midpoints for face ID: " << face->id() << std::endl;
            continue;
        }

        MeshLib::Vertex* m0 = edgeMidpoints[e0];
        MeshLib::Vertex* m1 = edgeMidpoints[e1];
        MeshLib::Vertex* m2 = edgeMidpoints[e2];

        if (vertexToID.find(v0) == vertexToID.end() ||
            vertexToID.find(v1) == vertexToID.end() ||
            vertexToID.find(v2) == vertexToID.end() ||
            vertexToID.find(m0) == vertexToID.end() ||
            vertexToID.find(m1) == vertexToID.end() ||
            vertexToID.find(m2) == vertexToID.end()) {
            std::cerr << "Error: Missing vertex IDs for face ID: " << face->id() << std::endl;
            continue;
		}

        int v0_id = vertexToID[v0];
        int v1_id = vertexToID[v1];
        int v2_id = vertexToID[v2];
        int m0_id = vertexToID[m0];
        int m1_id = vertexToID[m1];
        int m2_id = vertexToID[m2];

        int vid[3];

        // First face
        vid[0] = v0_id; vid[1] = m0_id; vid[2] = m2_id;
        Face* newFace1 = delegate.createFace(&newMesh, vid, nextFaceID++);
        if (newFace1 == nullptr) {
            std::cerr << "Failed to create face 1 for original face ID: " << face->id() << std::endl;
            continue;
        }
        faces_created++;

        // Second face
        vid[0] = v1_id; vid[1] = m1_id; vid[2] = m0_id;
        Face* newFace2 = delegate.createFace(&newMesh, vid, nextFaceID++);
        if (newFace2 == nullptr) {
            std::cerr << "Failed to create face 2 for original face ID: " << face->id() << std::endl;
            continue;
        }
        faces_created++;

        // Third face
        vid[0] = v2_id; vid[1] = m2_id; vid[2] = m1_id;
        Face* newFace3 = delegate.createFace(&newMesh, vid, nextFaceID++);
        if (newFace3 == nullptr) {
            std::cerr << "Failed to create face 3 for original face ID: " << face->id() << std::endl;
            continue;
        }
        faces_created++;

        // Central face
        vid[0] = m0_id; vid[1] = m1_id; vid[2] = m2_id;
        Face* newFace4 = delegate.createFace(&newMesh, vid, nextFaceID++);
        if (newFace4 == nullptr) {
            std::cerr << "Failed to create central face for original face ID: " << face->id() << std::endl;
            continue;
        }
        faces_created++;
    }

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

	printf("Write out the resultant obj file\n");
}
