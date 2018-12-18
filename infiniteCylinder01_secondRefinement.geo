//+
R_inner=1e-6;
R_outer=2e-6;
cellSize=.0125*(R_outer-R_inner);
Point(1) = {-R_outer, 0, 0, cellSize};
//+
Point(2) = {-R_inner, 0, 0, cellSize};
//+
Point(3) = {0, R_inner, 0, cellSize};
//+
Point(4) = {R_inner, 0, 0, cellSize};
//+
Point(5) = {R_outer, 0, 0, cellSize};
//+
Point(6) = {0, R_outer, 0, cellSize};
//+
Point(7) = {0, 0, 0, cellSize};
//+
Circle(1) = {4, 7, 2};
//+
Circle(2) = {5, 7, 1};
//+
Line(3) = {1, 2};
//+
Line(4) = {4, 5};
//+
Line Loop(1) = {3,-1,4,2};
//+
Plane Surface(1) = {1};
// at this point you have a uniform mesh of triangles
Recombine Surface{1};
// now you have a paved-looking collection of quads
// can also try the experimental 'delauney for quads' algorithm:
Mesh.Algorithm=8;
// Mesh.Smoothing=100;
