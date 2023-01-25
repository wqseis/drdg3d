// Gmsh project created on Sat Jul 23 16:52:46 2022
SetFactory("OpenCASCADE");

point_size = 1/Sqrt(3)*2;

Point(1) = {0, -20,   0, point_size};
Point(2) = {0,  20,   0, point_size};
Point(3) = {0,  20, -20, point_size};
Point(4) = {0, -20, -20, point_size};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
//Box(1) = {-50, -50, -50, 100, 100, 50};
Box(1) = {-100, -100, -100, 200, 200, 100};

// set size at domain boundary
MeshSize {5, 7, 11, 9, 10, 12, 8, 6} = 20/Sqrt(3)*2;

BooleanFragments{ Volume{1}; Delete; }{ Surface{1}; Delete; }

Mesh.Algorithm = 6; // Frontal-Delaunay for 2D meshes
//
