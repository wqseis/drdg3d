// Gmsh project created on Wed Jun  8 01:34:59 2022
SetFactory("OpenCASCADE");
//+
Box(1) = {-40, -50, -30, 80, 100, 60};
//+
MeshSize {3, 7, 5, 6, 8, 4, 2, 1} = 5;
//+
Point(9)  = {0, -15,  7.5, 0.2};
//+
Point(10) = {0,  15,  7.5, 0.2};
//+
Point(11) = {0,  15, -7.5, 0.2};
//+
Point(12) = {0, -15, -7.5, 0.2};
//+
Line(13) = {10, 9};
//+
Line(14) = {9, 12};
//+
Line(15) = {12, 11};
//+
Line(16) = {11, 10};
//+
Curve Loop(7) = {13, 14, 15, 16};
//+
Plane Surface(7) = {7};

Surface{7} In Volume{1};
//Line{13} In Surface{6};

Mesh.Algorithm = 6; // Frontal-Delaunay for 2D meshes
//

//

//
//+
