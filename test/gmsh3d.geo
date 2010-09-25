// Gmsh project created on Wed Aug  4 14:08:42 2010
cl = 0.1;
ccl = 0.5 * cl;
Point(1) = {0, 0, 0, cl};
Point(2) = {1, 0, 0, cl};
Point(3) = {1, 1, 0, cl};
Point(4) = {0, 1, 0, cl};
Point(5) = {0.2, 0, 0, ccl};
Point(6) = {0.7, 1, 0, ccl};
Point(7) = {0.65, 0.25, 0, ccl};
Point(8) = {0.4, 0.4, 0, ccl};
Point(9) = {0.4, 0.5, 0, ccl};
Point(10) = {0.6, 0.7, 0, ccl};
Point(11) = {0.8, 0.75, 0, ccl};
Point(12) = {0.25, 0.1, 0, ccl};

Line(1) = {1, 5};
Line(2) = {5, 2};
Line(3) = {2, 3};
Line(4) = {3, 6};
Line(5) = {6, 4};
Line(6) = {4, 1};

BSpline(9) = {5, 12, 7, 8, 9, 10, 11, 6};
Line Loop(10) = {2, 3, 4, -9};
Plane Surface(11) = {10};
Line Loop(12) = {6, 1, 9, 5};
Plane Surface(13) = {12};
Physical Surface(0) = {11};
Physical Surface(1) = {13};
