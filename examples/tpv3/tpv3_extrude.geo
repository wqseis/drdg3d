// Gmsh project created on Sat Jul 23 16:52:46 2022
SetFactory("OpenCASCADE");
//+
s = 1;
Point(1) = {0, -15, -7.5, s};
Point(2) = {0, -15, 7.5, s};
Point(3) = {0, 15, 7.5, s};
Point(4) = {0, 15, -7.5, s};
Line(1) = {3, 4};
Line(2) = {4, 1};
Line(3) = {1, 2};
Line(4) = {2, 3};
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
//+
Extrude {s, 0, 0} {
  Surface{1}; Layers {1}; 
}
//+
Extrude {-s, 0, 0} {
  Surface{1}; Layers {1}; 
}
//+
Box(3) = {-40, -50, -30, 80, 100, 60};
//+
BooleanFragments{ Volume{1}; Volume{2}; Delete; }{ Volume{3}; Delete; }
//+
MeshSize {13, 15, 19, 20, 17, 14, 18, 16} = 5;

Mesh.Algorithm = 6; // Frontal-Delaunay for 2D meshes
//
