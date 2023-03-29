// Gmsh project created on Sat Jul 23 16:52:46 2022
SetFactory("OpenCASCADE");

lc_fault = 1;
lc = 20;
lc_DistMin = 5*lc_fault;

Point(1) = {0, -15, -7.5, lc_fault};
Point(2) = {0,  15, -7.5, lc_fault};
Point(3) = {0,  15,  7.5, lc_fault};
Point(4) = {0, -15,  7.5, lc_fault};
Point(5) = {0, -1.5, -1.5, lc_fault};
Point(6) = {0,  1.5, -1.5, lc_fault};
Point(7) = {0,  1.5,  1.5, lc_fault};
Point(8) = {0, -1.5,  1.5, lc_fault};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};
Curve Loop(1) = {1, 2, 3, 4};
Curve Loop(2) = {5, 6, 7, 8};
Plane Surface(1) = {1};
Plane Surface(2) = {2};

BooleanFragments{ Surface{1}; Delete; }{ Surface{2}; Delete;  }

// //Box(1) = {-100, -150, -100, 200, 300, 200};
Box(1) = {-100, -120, -100, 200, 240, 200};
//  
// set size at domain boundary
//MeshSize {5, 7, 11, 9, 10, 12, 8, 6} = lc;
MeshSize {9, 10, 11, 12, 13, 14, 15, 16} = lc;

BooleanFragments{ Volume{1}; Delete; }{ Surface{2,3}; Delete; }


// Sizing:

// Field[1] = Distance;
// //Field[1].FacesList = {1};
// Field[1].FacesList = {2,3};
// 
// // Matheval field
// Field[2] = MathEval;
// //Field[2].F = Sprintf("0.02*F1 + 0.00001*F1^2 + %g", lc_fault);
// //Field[2].F = Sprintf("0.02*F1 +(F1/2e3)^2 + %g", lc_fault);
// //Field[2].F = Sprintf("0.02*F1 +(F1/0.5e3)^2 + %g", lc_fault);
// //Field[2].F = Sprintf("0.00*F1 +((F1-5e3)/0.1e3)^2 + %g", lc_fault);
// Field[2].F = Sprintf("0.0000*F1 +((F1-0e3)/0.1)^2 + %g", lc_fault);
// 
// //// Equivalent of propagation size on element
// Field[3] = Threshold;
// Field[3].IField = 1;
// Field[3].LcMin = lc_fault;
// Field[3].LcMax = lc;
// Field[3].DistMin = lc_DistMin;
// Field[3].DistMax = 125*lc_fault+0.001;
// 
// //
// Field[4] = Min;
// Field[4].FieldsList = {2,3};
// //Field[4].FieldsList = {2};
// Background Field = 4;

//
//
//
//+
//+
//+
