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
Point(5) = {0, -1.5, -6, point_size};
Point(6) = {0,  1.5, -6, point_size};
Point(7) = {0,  1.5, -9, point_size};
Point(8) = {0, -1.5, -9, point_size};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};
Curve Loop(2) = {5, 6, 7, 8};
Plane Surface(2) = {2};
Point(9)  = {0, -9, -6, point_size};
Point(10) = {0, -6, -6, point_size};
Point(11) = {0, -6, -9, point_size};
Point(12) = {0, -9, -9, point_size};
Line(9)  = {9, 10};
Line(10) = {10, 11};
Line(11) = {11, 12};
Line(12) = {12, 9};
Curve Loop(3) = {9, 10, 11, 12};
Plane Surface(3) = {3};
Point(13) = {0, 6, -6, point_size};
Point(14) = {0, 9, -6, point_size};
Point(15) = {0, 9, -9, point_size};
Point(16) = {0, 6, -9, point_size};
Line(13) = {13, 14};
Line(14) = {14, 15};
Line(15) = {15, 16};
Line(16) = {16, 13};
Curve Loop(4) = {13, 14, 15, 16};
Plane Surface(4) = {4};
BooleanFragments{ Surface{1}; Delete; }{ Surface{2,3,4}; Delete; }
Box(1) = {-50, -50, -50, 100, 100, 50};

// set size at domain boundary
MeshSize {17, 18, 19, 20, 21, 22, 23, 24} = 10;

BooleanFragments{ Volume{1}; Delete; }{ Surface{2,3,4,5}; Delete; }

Mesh.Algorithm = 6; // Frontal-Delaunay for 2D meshes
//
//+
