Merge "test.geo";
// Merge "test.msh";

Field[1] = Cylinder;
Field[1].Radius = 100e-6;
Field[1].VIn = 3e-6;
Field[1].VOut = 1e-3;
Field[2] = Cylinder;
Field[2].Radius = 500e-6;
Field[2].VIn = 10e-6;
Field[2].VOut = 30e-6;
Field[10] = Min;
Field[10].FieldsList = {1, 2};
Background Field = 10;

Mesh.CharacteristicLengthMin = 2e-6;
Mesh.CharacteristicLengthMax = 500e-6;
// Mesh.CharacteristicLengthExtendFromBoundary = 0;
// Mesh.CharacteristicLengthFactor = 1;
// Mesh.CharacteristicLengthFromCurvature = 0;
// Mesh.CharacteristicLengthFromPoints = 0;

Coherence Mesh;
// RefineMesh;
// CreateTopology;
Mesh.Algorithm = 6;
Mesh.RemeshAlgorithm = 1;
Mesh.RemeshParametrization = 0;
Mesh.AngleSmoothNormals = 1;
Mesh.RecombinationAlgorithm = 0;
Mesh.OptimizeNetgen = 1;
Mesh.ColorCarousel = 2; // by phys entity
// Mesh.RandomFactor = 1e-9;
// Geometry.HideCompounds = 0;

ll[] = Line "*";
For i In {0 : #ll[]-1}
//Compound Line(newl) = ll[i];
EndFor
ss[] = Surface "*";
For i In {0 : #ss[]-1}
//Compound Surface(news) = ss[i];
EndFor
