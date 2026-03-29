% Octave validation script - runs RSP and saves numerical outputs for comparison
% Usage: octave-cli --no-gui run_octave_validation.m

clear all; close all; clear functions;

% Add all paths; octave_compat must be first to override Utils/readOBJ etc.
addpath(genpath(pwd));

% Install optim package if not already installed (needed for quadprog)
try
  pkg load optim;
catch
  disp('Installing optim package (needed for quadprog)...');
  pkg install -forge optim;
  pkg load optim;
end

% octave_compat must be added AFTER pkg load to override quadprog etc.
addpath('octave_compat');

path_save = 'Results/';
path_data = 'Mesh/';

if ~exist(path_save, 'dir')
  mkdir(path_save);
end

meshes = {'pig', 'B36', 'SquareMyles'};
configs = {
  struct('name','pig',         'ff','curvature', 'hardedge',false, 'boundary',true, 'seamless',true, 'energy','alignment'),
  struct('name','B36',         'ff','smooth',    'hardedge',true,  'boundary',true, 'seamless',true, 'energy','distortion'),
  struct('name','SquareMyles', 'ff','trivial',   'hardedge',false, 'boundary',true, 'seamless',true, 'energy','chebyshev')
};

for ci = 1:length(configs)
  cfg = configs{ci};
  clear functions;  % Clear persistent caches (sort_triangles) between meshes
  addpath(genpath(pwd)); pkg load optim; addpath('octave_compat');
  fprintf('\n=== Processing mesh: %s (ff=%s, energy=%s) ===\n', cfg.name, cfg.ff, cfg.energy);

  mesh_name = cfg.name;
  frame_field_type = cfg.ff;
  ifhardedge = cfg.hardedge;
  ifboundary = cfg.boundary;
  ifseamless_const = cfg.seamless;
  energy_type = cfg.energy;

  % Energy weights
  weight = struct();
  if strcmp(energy_type, 'distortion')
    weight.w_conf_ar = 0.5;
  end
  weight.w_gradv = 1e-2;

  %% Load mesh
  [X,T] = readOBJ([path_data, mesh_name, '.obj']);

  % Rescale: area equals one
  area_tot = sum(sqrt(sum(cross(X(T(:,1),:) - X(T(:,2),:), X(T(:,1),:) - X(T(:,3),:),2).^2,2)))/2;
  X = X/sqrt(area_tot);

  % Preprocess geometry
  Src = MeshInfo(X, T);
  dec = dec_tri(Src);
  [param,Src,dec] = preprocess_ortho_param(Src, dec, ifboundary, ifhardedge, 40);

  fprintf('  Vertices: %d, Faces: %d, Edges: %d\n', Src.nv, Src.nf, Src.ne);

  %% Compute initial cross field
  if strcmp(frame_field_type, 'curvature')
    [omega,ang,sing,kappa,Curv] = compute_curvature_cross_field(Src, param, dec, 30, 1e-1);
    weight.aspect_ratio = ((abs(kappa(:,1)) + 1e-5)./(abs(kappa(:,2)) + 1e-5));
    t = exp(5);
    weight.aspect_ratio = max(min(weight.aspect_ratio, t), 1/t);
    weight.ang_dir = ang;
    if strcmp(energy_type, 'alignment')
      weight.w_ang   = 1;
      weight.w_ratio = 1;
    end
  elseif strcmp(frame_field_type, 'smooth')
    [omega,ang,sing] = compute_face_cross_field(Src, param, dec, 10);
  elseif strcmp(frame_field_type, 'trivial')
    sing = zeros(Src.nv,1);
    sing(param.idx_bound) = round(2*param.K(param.idx_bound)/pi)/4;
    om_cycle = param.Icycle*param.para_trans;
    om_cycle = om_cycle - 2*pi*round(4*om_cycle/(2*pi))/4;
    om_link = param.Ilink*param.para_trans;
    om_link = om_link - 2*pi*round(4*om_link/(2*pi))/4;
    [omega,ang,sing] = trivial_connection(Src, param, dec, ifboundary, ifhardedge, sing);
  end

  fprintf('  Singularities: %d positive, %d negative\n', sum(sing > 1/8), sum(sing < -1/8));

  %% Compute cross field jumps and build reduction matrix
  [Edge_jump,v2t,base_tri] = reduce_corner_var_2d(Src);
  [k21,Reduction] = reduction_from_ff2d(Src, param, ang, omega, Edge_jump, v2t);

  %% Optimize
  itmax = 200;
  ifplot = false;
  u = zeros(Src.nv,1);
  v = zeros(Src.nv,1);
  [u,v,ut,vt,om,angn,flag] = optimize_RSP(omega, ang, u, v, Src, param, dec, Reduction, energy_type, weight, ifplot, itmax);

  fprintf('  Optimization converged: %d\n', flag);

  %% Compute parametrization
  [SrcCut,dec_cut,Align,Rot] = mesh_to_disk_seamless(Src, param, angn, sing, k21, ifseamless_const, ifboundary, ifhardedge);
  [Xp,dX] = parametrization_from_scales(Src, SrcCut, dec_cut, param, angn, om, ut, vt, Align, Rot);

  %% Extract metrics
  disto = extract_scale_from_param(Xp, Src.X, Src.T, param, SrcCut.T, angn);
  curl_dX = sqrt(sum((dec_cut.d1p*dX).^2,2))./Src.area;

  num_flipped = sum(disto.detJ <= 0);
  fprintf('  Flipped triangles: %d / %d (%.1f%%)\n', num_flipped, Src.nf, 100*num_flipped/Src.nf);
  fprintf('  Mean log area distortion: %.6f\n', mean(log10(disto.area)));
  fprintf('  Mean abs log conformal: %.6f\n', mean(abs(log10(disto.conf))));
  fprintf('  Mean curl: %.6e\n', mean(curl_dX));
  fprintf('  UV range: [%.6f, %.6f] x [%.6f, %.6f]\n', min(Xp(:,1)), max(Xp(:,1)), min(Xp(:,2)), max(Xp(:,2)));

  %% Save numerical data
  out_prefix = [path_save, mesh_name, '_octave_'];
  save([out_prefix, 'Xp.txt'], 'Xp', '-ascii', '-double');
  save([out_prefix, 'ang.txt'], 'ang', '-ascii', '-double');
  save([out_prefix, 'sing.txt'], 'sing', '-ascii', '-double');
  save([out_prefix, 'u.txt'], 'u', '-ascii', '-double');
  save([out_prefix, 'v.txt'], 'v', '-ascii', '-double');
  save([out_prefix, 'omega.txt'], 'omega', '-ascii', '-double');

  fprintf('  Saved numerical outputs to %s*\n', out_prefix);
end

fprintf('\n=== All meshes processed successfully ===\n');
