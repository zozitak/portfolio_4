//+
SetFactory("OpenCASCADE");
Box(1) = {-0.5, -0.1, -0.1, 1, 0.2, 0.2};
//+
Physical Surface("anchor", 13) = {1};
//+
Physical Surface("force", 14) = {4};
//+
Transfinite Curve {11, 12, 10, 9} = 41 Using Progression 1;
//+
Transfinite Curve {3, 2, 1, 4, 7, 6, 5, 8} = 9 Using Progression 1;
//+
Transfinite Surface {1, 2, 3, 4, 5, 6};
//+
Recombine Surface {6, 1, 4, 5, 2, 3};
//+
Transfinite Volume{1} = {1, 2, 4, 3, 5, 6, 8, 7};
