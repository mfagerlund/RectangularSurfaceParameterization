% Minimal Octave test - step by step
clear all; close all;
addpath(genpath(pwd));
addpath('octave_compat');
pkg load optim;

fprintf('Step 1: Loading mesh...\n');
[X,T] = readOBJ('Mesh/B36.obj');
fprintf('  Loaded: %d vertices, %d faces\n', size(X,1), size(T,1));

fprintf('Step 2: Rescaling...\n');
area_tot = sum(sqrt(sum(cross(X(T(:,1),:) - X(T(:,2),:), X(T(:,1),:) - X(T(:,3),:),2).^2,2)))/2;
X = X/sqrt(area_tot);
fprintf('  Area: %.6f\n', area_tot);

fprintf('Step 3: MeshInfo...\n');
Src = MeshInfo(X, T);
fprintf('  nv=%d, nf=%d, ne=%d\n', Src.nv, Src.nf, Src.ne);

fprintf('Step 4: dec_tri...\n');
dec = dec_tri(Src);
fprintf('  Done\n');

fprintf('Step 5: preprocess_ortho_param...\n');
[param,Src,dec] = preprocess_ortho_param(Src, dec, true, false, 40);
fprintf('  Done. ide_fix=%d\n', length(param.ide_fix));

fprintf('Step 6: compute_face_cross_field (smooth)...\n');
[omega,ang,sing] = compute_face_cross_field(Src, param, dec, 10);
fprintf('  Done. Singularities: %d pos, %d neg\n', sum(sing > 1/8), sum(sing < -1/8));

fprintf('All steps passed!\n');
