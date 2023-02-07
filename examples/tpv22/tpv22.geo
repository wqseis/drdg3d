//+
SetFactory("OpenCASCADE");

faultsize = 1.6/1.0/Sqrt(3)*2.0;
//faultsize = 1.6/8.0/Sqrt(3)*2.0;

//+
Point(1) = {0, -25, 0, faultsize};
Point(2) = {0, 5, 0, faultsize};
Point(3) = {0, 5, -20, faultsize};
Point(4) = {0, -25, -20, faultsize};
//+
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
//+
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
//+
Point(5) = {1.6, -5, 0, faultsize};
Point(6) = {1.6, 25, 0, faultsize};
Point(7) = {1.6, 25, -20, faultsize};
Point(8) = {1.6, -5, -20, faultsize};
//+
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};
//+
Curve Loop(2) = {5, 6, 7, 8};
Plane Surface(2) = {2};
//+
////Box(1) = {-80, -70, -80, 160, 200, 80};
//Box(1) = {-150, -150, -200, 300, 300, 200};
Box(1) = {-100, -100, -100, 200, 200, 100};
////+
BooleanFragments{ Volume{1}; Delete; }{ Surface{1}; Delete; }
BooleanFragments{ Volume{1}; Delete; }{ Surface{2}; Delete; }
////+
////MeshSize {:} = 20;
//
//MeshSize {:} = 40/Sqrt(3.0)*2.0;
//
Mesh.Algorithm = 6; // Frontal-Delaunay for 2D meshes
//+
MeshSize {9, 10, 11, 12, 13, 14, 15, 18} = 10/Sqrt(3)*2;
