// Gmsh project created on Mon Apr 11 21:52:13 2022

SetFactory("OpenCASCADE");

Merge "surf_400m.step";

//+
Rotate {{0, 1, 0}, {0, 0, 0}, Pi/2} {
  Surface{1};
}
//+
Rotate {{1, 0, 0}, {0, 0, 0}, Pi/2} {
  Surface{1};
}

////Box(1) = {-40, -40, -40, 80, 40, 80};
Box(1) = {-100, -100, -100, 200, 200, 100};
//
Surface{1} In Volume{1};
Curve{3} In Surface{7};

//
MeshSize {3, 4, 2, 1} = 1/Sqrt(3)*2/1;
//
MeshSize {5,6,7,8,9,10,11,12} = 20/Sqrt(3)*2;
//
Mesh.Algorithm = 6; // Frontal-Delaunay for 2D meshes
////


//+
Show "*";
