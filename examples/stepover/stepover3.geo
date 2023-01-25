//+
SetFactory("OpenCASCADE");

faultsize = 1.0/1.0/Sqrt(3)*2.0;

//+
Point(1) = {0, 0, -1, faultsize};
Point(2) = {0, 40, -1, faultsize};
Point(3) = {0, 40, -16, faultsize};
Point(4) = {0, 0, -16, faultsize};
//+
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
//+
Curve Loop(1) = {2, 3, 4, 1};
Plane Surface(1) = {1};
//+
Point(5) = {2, 32, -1, faultsize};
Point(6) = {2, 62, -1, faultsize};
Point(7) = {2, 62, -16, faultsize};
Point(8) = {2, 32, -16, faultsize};
//+
Line(5) = {6, 5};
Line(6) = {5, 8};
Line(7) = {8, 7};
Line(8) = {7, 6};
//+
Curve Loop(2) = {8, 5, 6, 7};
Plane Surface(2) = {2};
//+
Point(9)  = {1, 32, -1, faultsize};
Point(10) = {1, 40, -1, faultsize};
Point(11) = {1, 40, -16, faultsize};
Point(12) = {1, 32, -16, faultsize};
//+
Line(9) = {10, 11};
Line(10) = {11, 12};
Line(11) = {12, 9};
Line(12) = {9, 10};
//+
Curve Loop(3) = {9, 10, 11, 12};
Plane Surface(3) = {3};
//+
//Box(1) = {-80, -70, -80, 160, 200, 80};
Box(1) = {-150, -120, -200, 300, 300, 200};
//+
BooleanFragments{ Volume{1}; Delete; }{ Surface{1}; Delete; }
BooleanFragments{ Volume{1}; Delete; }{ Surface{2}; Delete; }
BooleanFragments{ Volume{1}; Delete; }{ Surface{3}; Delete; }
//+
//MeshSize {15, 19, 18, 13, 14, 17, 24, 16} = 20;

Point(101) = {0, 0, 0, faultsize};
Point(102) = {0, 40, 0, faultsize};
Point(103) = {1, 32, 0, faultsize};
Point(104) = {1, 40, 0, faultsize};
Point(105) = {2, 32, 0, faultsize};
Point(106) = {2, 62, 0, faultsize};
//+
Line(33) = {101, 102};
Line(34) = {103, 104};
Line(35) = {105, 106};
Line{33,34,35} In Surface{6};
//+
MeshSize {15, 16, 19, 20, 17, 18, 13, 14} = 40/Sqrt(3.0)*2.0;


Mesh.Algorithm = 6; // Frontal-Delaunay for 2D meshes
