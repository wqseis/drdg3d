// Gmsh project created on Sat Jul 23 16:52:46 2022
SetFactory("OpenCASCADE");

lc_fault = 0.2;
lc = 10;
lc_DistMin = 2*lc_fault;

Point(1) = {0, -20,   0, lc_fault};
Point(2) = {0,  20,   0, lc_fault};
Point(3) = {0,  20, -20, lc_fault};
Point(4) = {0, -20, -20, lc_fault};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
//Box(1) = {-50, -50, -50, 100, 100, 50};
Box(1) = {-100, -100, -100, 200, 200, 100};

// set size at domain boundary
MeshSize {7, 5, 6, 8, 11, 12, 10, 9} = lc;

BooleanFragments{ Volume{1}; Delete; }{ Surface{1}; Delete; }

// Sizing:

Field[1] = Distance;
Field[1].FacesList = {1};

// Matheval field
Field[2] = MathEval;
//Field[2].F = Sprintf("0.02*F1 + 0.00001*F1^2 + %g", lc_fault);
//Field[2].F = Sprintf("0.02*F1 +(F1/2e3)^2 + %g", lc_fault);
//Field[2].F = Sprintf("0.02*F1 +(F1/0.5e3)^2 + %g", lc_fault);
//Field[2].F = Sprintf("0.00*F1 +((F1-5e3)/0.1e3)^2 + %g", lc_fault);
Field[2].F = Sprintf("0.0000*F1 +((F1-0e3)/5)^2 + %g", lc_fault);

//// Equivalent of propagation size on element
Field[3] = Threshold;
Field[3].IField = 1;
Field[3].LcMin = lc_fault;
Field[3].LcMax = lc;
Field[3].DistMin = lc_DistMin;
Field[3].DistMax = 50+0.001;

//
Field[4] = Min;
Field[4].FieldsList = {2,3};
//Field[4].FieldsList = {2};
Background Field = 4;

Mesh.Algorithm = 6; // Frontal-Delaunay for 2D meshes
//
//+
//+
