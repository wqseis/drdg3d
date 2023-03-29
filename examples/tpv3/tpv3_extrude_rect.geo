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

Point(5) = {0, -1.5, -1.5, s};
Point(6) = {0, -1.5, +1.5, s};
Point(7) = {0,  1.5, +1.5, s};
Point(8) = {0,  1.5, -1.5, s};
Line(5) = {7, 8};
Line(6) = {8, 5};
Line(7) = {5, 6};
Line(8) = {6, 7};
Curve Loop(2) = {5, 6, 7, 8};
Plane Surface(2) = {2};


Mesh.Algorithm = 6; // Frontal-Delaunay for 2D meshes
//
//+
BooleanFragments{ Surface{1}; Delete; }{Surface{2}; Delete;  }
//+
Extrude {s, 0, 0} {
  Surface{2, 3}; Layers {1};
}
Extrude {-s, 0, 0} {
  Surface{2, 3}; Layers {1};
}

Box(5) = {-40, -50, -30, 80, 100, 60};
BooleanFragments{ Volume{1, 2, 3, 4}; Delete; }{ Volume{3}; Delete; }
//+
//+
MeshSize {27, 25, 29, 31, 26, 30, 32, 28} = 5;
