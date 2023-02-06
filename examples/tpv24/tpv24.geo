// Gmsh project created on Sun Aug  7 09:23:16 2022
SetFactory("OpenCASCADE");

s = 1.0/2/Sqrt(3.0)*2.0;
s = 0.5;

Point(1) = {0, -16, 0, s};
Point(2) = {0, -16, -15, s};
Point(3) = {0,  12, 0, s};
Point(4) = {0,  12, -15, s};

ang=30/180*Pi;
Point(5) = {0, 0, -15, s};
Point(6) = {0, 0, 0, s};
Point(7) = {12*Sin(ang), 12*Cos(ang), -15, s};
Point(8) = {12*Sin(ang), 12*Cos(ang), 0, s};

Line(13) = {1, 6};
//+
Line(14) = {6, 3};
//+
Line(15) = {3, 4};
//+
Line(16) = {4, 5};
//+
Line(17) = {5, 2};
//+
Line(18) = {2, 1};
//+
Line(19) = {6, 5};
//+
Line(20) = {5, 7};
//+
Line(21) = {7, 8};
//+
Line(22) = {8, 6};
//+
Curve Loop(1) = {13, 14, 15, 16, 17, 18};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {22, 19, 20, 21};
//+
Plane Surface(2) = {2};

//+
//Box(1) = {-80, -80, -80, 160, 160, 80};
Box(1) = {-100, -100, -120, 200, 200, 120};


//+
BooleanFragments{ Volume{1}; Delete; }{ Surface{1}; Surface{2}; Delete; }
//+
MeshSize {15, 13, 14, 16, 19, 20, 17, 18} = 20/Sqrt(3)*2;
