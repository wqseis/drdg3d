// Gmsh project created on Sat Jul 23 16:52:46 2022
SetFactory("OpenCASCADE");

lc_fault = 1;
lc = 50;
lc_DistMin = 1*lc_fault;

// for benchmark
//lc_fault = 0.2;
//lc = 10;
//lc_DistMin = 5*lc_fault;


Point(1) = {0, -15,   0, lc_fault};
Point(2) = {0,  15,   0, lc_fault};
Point(3) = {0,  15, -15, lc_fault};
Point(4) = {0, -15, -15, lc_fault};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Point(5) = {0, -1.5, -6, lc_fault};
Point(6) = {0,  1.5, -6, lc_fault};
Point(7) = {0,  1.5, -9, lc_fault};
Point(8) = {0, -1.5, -9, lc_fault};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};
Curve Loop(2) = {5, 6, 7, 8};
Plane Surface(2) = {2};
Point(9)  = {0, -9, -6, lc_fault};
Point(10) = {0, -6, -6, lc_fault};
Point(11) = {0, -6, -9, lc_fault};
Point(12) = {0, -9, -9, lc_fault};
Line(9)  = {9, 10};
Line(10) = {10, 11};
Line(11) = {11, 12};
Line(12) = {12, 9};
Curve Loop(3) = {9, 10, 11, 12};
Plane Surface(3) = {3};
Point(13) = {0, 6, -6, lc_fault};
Point(14) = {0, 9, -6, lc_fault};
Point(15) = {0, 9, -9, lc_fault};
Point(16) = {0, 6, -9, lc_fault};
Line(13) = {13, 14};
Line(14) = {14, 15};
Line(15) = {15, 16};
Line(16) = {16, 13};
Curve Loop(4) = {13, 14, 15, 16};
Plane Surface(4) = {4};
BooleanFragments{ Surface{1}; Delete; }{ Surface{2,3,4}; Delete; }
//Box(1) = {-50, -50, -50, 100, 100, 50};
Box(1) = {-100, -100, -100, 200, 200, 100};

// set size at domain boundary
MeshSize {17, 18, 19, 20, 21, 22, 23, 24} = lc;

BooleanFragments{ Volume{1}; Delete; }{ Surface{2,3,4,5}; Delete; }

// Sizing:

Field[1] = Distance;
Field[1].FacesList = {2,3,4,5};

// Matheval field
Field[2] = MathEval;
//Field[2].F = Sprintf("0.02*F1 + 0.00001*F1^2 + %g", lc_fault);
//Field[2].F = Sprintf("0.02*F1 +(F1/2e3)^2 + %g", lc_fault);
//Field[2].F = Sprintf("0.02*F1 +(F1/0.5e3)^2 + %g", lc_fault);
//Field[2].F = Sprintf("0.00*F1 +((F1-5e3)/0.1e3)^2 + %g", lc_fault);
Field[2].F = Sprintf("0.0000*F1 +((F1-0e3)/0.1)^2 + %g", lc_fault);
// for benchmark
//Field[2].F = Sprintf("0.0000*F1 +((F1-0e3)/5)^2 + %g", lc_fault);

//// Equivalent of propagation size on element
Field[3] = Threshold;
Field[3].IField = 1;
Field[3].LcMin = lc_fault;
Field[3].LcMax = lc;
Field[3].DistMin = lc_DistMin;
Field[3].DistMax = 125*lc_fault+0.001;

//
Field[4] = Min;
Field[4].FieldsList = {2,3};
Background Field = 4;


Mesh.Algorithm = 6; // Frontal-Delaunay for 2D meshes
//
//+
