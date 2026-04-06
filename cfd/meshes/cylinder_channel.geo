Lx = 4.0;
Ly = 2.0;
cx = 1.5;
cy = 1.0;
R  = 0.25;

h_wall = 0.01;
h_far  = 0.1;

Point(1) = {0,  0,  0, h_far};
Point(2) = {Lx, 0,  0, h_far};
Point(3) = {Lx, Ly, 0, h_far};
Point(4) = {0,  Ly, 0, h_far};

Point(5) = {cx,   cy,   0, h_wall};
Point(6) = {cx+R, cy,   0, h_wall};
Point(7) = {cx,   cy+R, 0, h_wall};
Point(8) = {cx-R, cy,   0, h_wall};
Point(9) = {cx,   cy-R, 0, h_wall};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Circle(5) = {6, 5, 7};
Circle(6) = {7, 5, 8};
Circle(7) = {8, 5, 9};
Circle(8) = {9, 5, 6};

Curve Loop(1) = {1, 2, 3, 4};
Curve Loop(2) = {5, 6, 7, 8};
Plane Surface(1) = {1, 2};

Physical Curve("inflow", 10)   = {4};
Physical Curve("outflow", 11)  = {2};
Physical Curve("symmetry", 12) = {1, 3};
Physical Curve("wall", 13)     = {5, 6, 7, 8};
Physical Surface("fluid", 1)  = {1};

Field[1] = Distance;
Field[1].CurvesList = {5, 6, 7, 8};
Field[1].Sampling = 200;

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = h_wall;
Field[2].SizeMax = h_far;
Field[2].DistMin = 0.0;
Field[2].DistMax = 1.0;

Background Field = 2;

Mesh.Algorithm = 5;
