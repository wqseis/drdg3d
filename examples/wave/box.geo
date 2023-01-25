// Gmsh project created on Sat Jul 23 16:52:46 2022
SetFactory("OpenCASCADE");

Box(1) = {-50, -50, -50, 100, 100, 50};
 
// set size at domain boundary
MeshSize {:} = 5/Sqrt(3)*2;


//
//
