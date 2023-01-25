// Gmsh project created on Mon Jul 18 15:46:01 2022
SetFactory("OpenCASCADE");
//+
Box(1) = {-80, -80, -80, 160, 160, 80};
Rectangle(7) = {-15, -15, 0, 30, 30, 0};

DipAngle = 15;
Rotate {{0, 1, 0}, {0, 0, 0}, DipAngle/180.0*Pi} {
  Surface{7};
}

BooleanFragments{ Volume{1}; Delete; }{ Surface{7}; Delete; }
Recursive Delete {
  Surface{1};
}

MeshSize{:} = 20;
Characteristic Length {3, 5, 6, 4} = 0.75;
//MeshSize {9, 7, 12, 13, 14, 11, 10, 8} = 30;

Mesh.Algorithm = 6; // Frontal-Delaunay for 2D meshes
