// Gmsh project created on Sat Jul 29 10:20:43 2023
SetFactory("OpenCASCADE");
Mesh.MshFileVersion = 2.2;
//+
Cylinder(1) = {0, 0, 0, 0, 0, 30, 2, 2*Pi};
//+
Physical Surface("inlet", 1) = {3};
//+
Physical Surface("outlet", 2) = {2};
//+
Physical Surface("wall", 3) = {1};
//+
Physical Volume("fluid", 1) = {1};
//+
Field[1] = Constant;
//+
Field[1].CurvesList = {1, 2, 3};
Field[1].IncludeBoundary = 0;
//+
Field[1].VIn = 0.15;
//+
Field[1].VOut = 0.15;
//+
Field[1].VolumesList = {1};
//+
Background Field = 1;

Mesh 3;
