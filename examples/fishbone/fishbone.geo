// Gmsh project created on Sun Aug  7 09:23:16 2022
SetFactory("OpenCASCADE");

s=1.0/Sqrt(3.0)*2.0;

Point(1) = {0, -20, 0, s};
Point(2) = {0, -20, -15, s};
Point(3) = {0,  20, 0, s};
Point(4) = {0,  20, -15, s};
 

Point(5) = {0, 0, -15, s};
Point(6) = {0, 0, 0, s};
Point(7) = {5, 5, -15, s};
Point(8) = {5, 5, 0, s};
Point(9)  = {0, 5, -15, s};
Point(10) = {0, 5, 0, s};
Point(11) = {5, 10, -15, s};
Point(12) = {5, 10, 0, s};
Point(13) = {0, 10, -15, s};
Point(14) = {0, 10, 0, s};
Point(15) = {5, 15, -15, s};
Point(16) = {5, 15, 0, s};

//+
Line(1) = {1, 2};
//+
Line(2) = {2, 5};
//+
Line(3) = {5, 9};
//+
Line(4) = {9, 13};
//+
Line(5) = {13, 4};
//+
Line(6) = {4, 3};
//+
Line(7) = {3, 14};
//+
Line(8) = {14, 10};
//+
Line(9) = {10, 6};
//+
Line(10) = {6, 1};
//+
Line(11) = {6, 5};
//+
Line(12) = {10, 9};
//+
Line(13) = {14, 13};
//+
Line(14) = {8, 7};
//+
Line(15) = {5, 7};
//+
Line(16) = {6, 8};
//+
Line(17) = {10, 12};
//+
Line(18) = {12, 11};
//+
Line(19) = {11, 9};
//+
Line(20) = {14, 16};
//+
Line(21) = {16, 15};
//+
Line(22) = {15, 13};
//+
Curve Loop(1) = {10, 1, 2, -11};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {9, 11, 3, -12};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {8, 12, 4, -13};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {13, 5, 6, 7};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {20, 21, 22, -13};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {17, 18, 19, -12};
//+
Plane Surface(6) = {6};
//+
Curve Loop(7) = {16, 14, -15, -11};
//+
Plane Surface(7) = {7};
//+
Box(1) = {-80, -80, -80, 160, 160, 80};
//+
BooleanFragments{ Volume{1}; Delete; }{ Surface{1,2,3,4,5,6,7}; Delete; }
//BooleanFragments{ Volume{1}; Delete; }{ Surface{2}; Delete; }
//BooleanFragments{ Volume{1}; Delete; }{ Surface{3}; Delete; }
//BooleanFragments{ Volume{1}; Delete; }{ Surface{4}; Delete; }
//BooleanFragments{ Volume{1}; Delete; }{ Surface{13}; Delete; }
//BooleanFragments{ Volume{1}; Delete; }{ Surface{5}; Delete; }
//BooleanFragments{ Volume{1}; Delete; }{ Surface{16}; Delete; }
//+
MeshSize {19, 23, 20, 24, 22, 21, 17, 18} = 20/Sqrt(3.0)*2.0;
