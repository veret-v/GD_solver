
SetFactory("OpenCASCADE");

lc = 0.035;
lc_step = 0.018;

Point(1) = {0.0, 0.2, 0, lc_step};
Point(2) = {0.6, 0.2, 0, lc_step};
Point(3) = {0.6, 0.0, 0, lc_step};
Point(4) = {3.0, 0.0, 0, lc};
Point(5) = {3.0, 1.0, 0, lc};
Point(6) = {0.0, 1.0, 0, lc};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 1};
Curve Loop(1) = {1, 2, 3, 4, 5, 6};
Plane Surface(1) = {1};

Field[1] = Box;
Field[1].XMin = 0.45; Field[1].XMax = 1.8;
Field[1].YMin = 0.0;  Field[1].YMax = 0.7;
Field[1].VIn = 0.018;
Field[1].VOut = lc;
Background Field = 1;

Physical Curve("wall", 10) = {1, 2, 3, 5};
Physical Curve("outflow", 11) = {4};
Physical Curve("inflow", 12) = {6};
Physical Surface("fluid", 20) = {1};

Mesh.Algorithm = 5;
Mesh.CharacteristicLengthMin = 0.014;
Mesh.CharacteristicLengthMax = 0.045;
Mesh.Optimize = 1;
