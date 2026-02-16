// Gmsh project created on Wed Oct 19 10:54:38 2022
SetFactory("OpenCASCADE");
//+
Point(1) = {-0.5, -0.5, 0, 1.0};
//+
Point(2) = {-0.5, 0, 0, 1.0};
//+
Point(3) = {-0.5, 0.5, 0, 1.0};
//+
Point(4) = {0, -0.5, 0, 1.0};
//+
Point(5) = {0.5, -0.5, 0, 1.0};
//+
Point(6) = {0.5, 0, 0, 1.0};
//+
Point(7) = {0.5, 0.5, 0, 1.0};
//+
Point(8) = {0, 0.5, 0, 1.0};
//+
Point(9) = {0, 0, 0, 1.0};
//+
Line(1) = {1, 4};
//+
Line(2) = {4, 5};
//+
Line(3) = {5, 6};
//+
Line(4) = {7, 6};
//+
Line(5) = {8, 7};
//+
Line(6) = {3, 8};
//+
Line(7) = {2, 3};
//+
Line(8) = {1, 2};
//+
Line(9) = {2, 9};
//+
Line(10) = {9, 8};
//+
Line(11) = {4, 9};
//+
Line(12) = {9, 6};
//+
Curve Loop(1) = {7, 6, -10, -9};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {5, 4, -12, 10};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {11, 12, -3, -2};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {9, -11, -1, 8};
//+
Plane Surface(4) = {4};
//+
Physical Curve("Inflow", 1) = {6, 4, 2, 8};
//+
Physical Curve("Outflow", 2) = {7, 5, 3, 1};
//+
Physical Surface(14) = {1};
//+
Physical Surface(15) = {2};
//+
Physical Surface(16) = {3};
//+
Physical Surface(17) = {4};

Mesh.MshFileVersion = 2.2;
