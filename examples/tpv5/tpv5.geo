// Gmsh project created on Sat Jul 23 16:52:46 2022
SetFactory("OpenCASCADE");

point_size = 1;///Sqrt(3)*2;

Point(1) = {0, -15,   0, point_size};
Point(2) = {0,  15,   0, point_size};
Point(3) = {0,  15, -15, point_size};
Point(4) = {0, -15, -15, point_size};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Box(1) = {-50, -50, -50, 100, 100, 50};

// set size at domain boundary
MeshSize {5, 7, 11, 9, 10, 12, 8, 6} = 10;

BooleanFragments{ Volume{1}; Delete; }{ Surface{1}; Delete; }

Mesh.Algorithm = 6; // Frontal-Delaunay for 2D meshes
//
