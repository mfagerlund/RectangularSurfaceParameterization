clear all; close all;
addpath(genpath(pwd));
pkg load optim;
addpath('octave_compat');

mesh_name = 'pig';
fprintf('Loading %s...\n', mesh_name);
[X,T] = readOBJ(['Mesh/', mesh_name, '.obj']);
area_tot = sum(sqrt(sum(cross(X(T(:,1),:) - X(T(:,2),:), X(T(:,1),:) - X(T(:,3),:),2).^2,2)))/2;
X = X/sqrt(area_tot);
Src = MeshInfo(X, T);
dec = dec_tri(Src);
fprintf('  nv=%d, nf=%d, ne=%d\n', Src.nv, Src.nf, Src.ne);

fprintf('Preprocessing...\n');
[param,Src,dec] = preprocess_ortho_param(Src, dec, true, false, 40);

fprintf('Smooth cross field...\n');
[omega,ang,sing] = compute_face_cross_field(Src, param, dec, 10);
fprintf('  Singularities: %d pos, %d neg\n', sum(sing > 1/8), sum(sing < -1/8));

fprintf('Reduction...\n');
[Edge_jump,v2t,base_tri] = reduce_corner_var_2d(Src);
[k21,Reduction] = reduction_from_ff2d(Src, param, ang, omega, Edge_jump, v2t);

fprintf('Optimizing (distortion, w_conf_ar=0.5)...\n');
weight.w_conf_ar = 0.5;
weight.w_gradv = 1e-2;
u = zeros(Src.nv,1);
v = zeros(Src.nv,1);
[u,v,ut,vt,om,angn,flag] = optimize_RSP(omega, ang, u, v, Src, param, dec, Reduction, 'distortion', weight, false, 200);

fprintf('Computing parametrization...\n');
[SrcCut,dec_cut,Align,Rot] = mesh_to_disk_seamless(Src, param, angn, sing, k21, true, true, false);
[Xp,dX] = parametrization_from_scales(Src, SrcCut, dec_cut, param, angn, om, ut, vt, Align, Rot);

disto = extract_scale_from_param(Xp, Src.X, Src.T, param, SrcCut.T, angn);
num_flipped = sum(disto.detJ <= 0);
fprintf('  Flipped: %d / %d (%.1f%%)\n', num_flipped, Src.nf, 100*num_flipped/Src.nf);
fprintf('  UV range: [%.6f, %.6f] x [%.6f, %.6f]\n', min(Xp(:,1)), max(Xp(:,1)), min(Xp(:,2)), max(Xp(:,2)));

save('Results/pig_octave_Xp.txt', 'Xp', '-ascii', '-double');
save('Results/pig_octave_ang.txt', 'ang', '-ascii', '-double');
save('Results/pig_octave_sing.txt', 'sing', '-ascii', '-double');
save('Results/pig_octave_u.txt', 'u', '-ascii', '-double');
save('Results/pig_octave_v.txt', 'v', '-ascii', '-double');
save('Results/pig_octave_omega.txt', 'omega', '-ascii', '-double');
fprintf('DONE\n');
