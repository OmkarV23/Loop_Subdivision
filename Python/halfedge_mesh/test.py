''' read the file line by line to understanding how the half edge load 
an off file and to save the half edge mesh as an obj file.

author:
    zhangsihao yang
'''
print 'hello'

data_path = 'tests/data/cube.off'

import halfedge_mesh

# HalfedgeMesh
mesh = halfedge_mesh.HalfedgeMesh(data_path)

# Returns a list of Vertex type (in order of file)--similarly for halfedges,
# and facets
mesh.vertices

# The number of facets in the mesh
print(len(mesh.facets))

# Get the 10th halfedge
mesh.halfedges[10]

# Get the halfedge that starts at vertex 25 and ends at vertex 50
# mesh.get_halfedge(0, 1)

print(mesh.vertices)
for vertex in mesh.vertices:
    print(vertex.get_vertex())

print(mesh.facets)
print('--------------------')
for face in mesh.facets:
    print(face.a, face.b, face.c)

# to save the halfedge mesh you will need to following function.
def save_halfmesh_as_obj(mesh, file_name):
    with open(file_name, 'w') as open_file:
        for vertex in mesh.vertices:
            lv = vertex.get_vertex()
            open_file.write("v {} {} {} \n".format(lv[0], lv[1], lv[2]))

        for face in mesh.facets:
            open_file.write("f {} {} {}\n".format(face.a+1, face.b+1, face.c+1))

save_halfmesh_as_obj(mesh, 'tmp.obj')

