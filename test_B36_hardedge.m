clear all; close all; clear functions;
addpath(genpath(pwd));
pkg load optim;
addpath('octave_compat');

mesh_name = 'B36';
fprintf('Loading %s...\n', mesh_name);
[X,T] = readOBJ(['Mesh/', mesh_name, '.obj']);
area_tot = sum(sqrt(sum(cross(X(T(:,1),:) - X(T(:,2),:), X(T(:,1),:) - X(T(:,3),:),2).^2,2)))/2;
X = X/sqrt(area_tot);
Src = MeshInfo(X, T);
dec = dec_tri(Src);

fprintf('Running preprocess with hardedge=true...\n');
try
  [param,Src,dec] = preprocess_ortho_param(Src, dec, true, true, 40);
  fprintf('SUCCESS: nv=%d, nf=%d, ne=%d\n', Src.nv, Src.nf, Src.ne);
catch e
  fprintf('ERROR: %s\n', e.message);
  % Debug: find which vertex causes the issue
  fprintf('Debugging sort_triangles for boundary vertices...\n');
end
