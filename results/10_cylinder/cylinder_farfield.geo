
SetFactory("OpenCASCADE");
R = 1.0;
Rfar = 15.0;
lc_body = 0.045;
lc_far = 0.75;

Point(1) = {0, 0, 0, lc_body};
Point(2) = { R, 0, 0, lc_body};
Point(3) = {0,  R, 0, lc_body};
Point(4) = {-R, 0, 0, lc_body};
Point(5) = {0, -R, 0, lc_body};
Point(6) = { Rfar, 0, 0, lc_far};
Point(7) = {0,  Rfar, 0, lc_far};
Point(8) = {-Rfar, 0, 0, lc_far};
Point(9) = {0, -Rfar, 0, lc_far};

Circle(1) = {2,1,3};
Circle(2) = {3,1,4};
Circle(3) = {4,1,5};
Circle(4) = {5,1,2};
Circle(5) = {6,1,7};
Circle(6) = {7,1,8};
Circle(7) = {8,1,9};
Circle(8) = {9,1,6};

Curve Loop(1) = {5,6,7,8};
Curve Loop(2) = {1,2,3,4};
Plane Surface(1) = {1,2};

Field[1] = Distance;
Field[1].CurvesList = {1,2,3,4};
Field[1].Sampling = 100;
Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = lc_body;
Field[2].SizeMax = lc_far;
Field[2].DistMin = 0.0;
Field[2].DistMax = 6.0;
Background Field = 2;

Physical Curve("cylinder", 10) = {1,2,3,4};
Physical Curve("farfield", 11) = {5,6,7,8};
Physical Surface("fluid", 20) = {1};

Mesh.Algorithm = 5;
Mesh.CharacteristicLengthMin = 0.035;
Mesh.CharacteristicLengthMax = 0.85;
Mesh.Optimize = 1;
