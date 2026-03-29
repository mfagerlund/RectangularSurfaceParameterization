function [V, F, UV, TF, N, NF, SI] = readOBJ(filename, varargin)
% READOBJ  Fast OBJ reader for Octave. Replaces the slow fscanf-based version.
  text = fileread(filename);
  lines = strsplit(text, '\n');

  verts = {};
  faces = {};
  uvs = {};
  tfaces = {};
  norms = {};
  nfaces = {};

  for i = 1:length(lines)
    ln = strtrim(lines{i});
    if isempty(ln) || ln(1) == '#'
      continue;
    end
    tokens = strsplit(ln);
    tag = tokens{1};

    if strcmp(tag, 'v') && length(tokens) >= 4
      verts{end+1} = [str2double(tokens{2}), str2double(tokens{3}), str2double(tokens{4})];
    elseif strcmp(tag, 'vt') && length(tokens) >= 3
      uvs{end+1} = [str2double(tokens{2}), str2double(tokens{3})];
    elseif strcmp(tag, 'vn') && length(tokens) >= 4
      norms{end+1} = [str2double(tokens{2}), str2double(tokens{3}), str2double(tokens{4})];
    elseif strcmp(tag, 'f')
      fv = [];
      ft = [];
      fn = [];
      for j = 2:length(tokens)
        parts = strsplit(tokens{j}, '/');
        fv(end+1) = str2double(parts{1});
        if length(parts) >= 2 && ~isempty(parts{2})
          ft(end+1) = str2double(parts{2});
        end
        if length(parts) >= 3 && ~isempty(parts{3})
          fn(end+1) = str2double(parts{3});
        end
      end
      % Triangulate if needed
      for j = 2:length(fv)-1
        faces{end+1} = [fv(1), fv(j), fv(j+1)];
        if ~isempty(ft)
          tfaces{end+1} = [ft(1), ft(j), ft(j+1)];
        end
        if ~isempty(fn)
          nfaces{end+1} = [fn(1), fn(j), fn(j+1)];
        end
      end
    end
  end

  V = vertcat(verts{:});
  if ~isempty(faces)
    F = vertcat(faces{:});
  else
    F = zeros(0, 3);
  end
  if ~isempty(uvs)
    UV = vertcat(uvs{:});
  else
    UV = zeros(0, 2);
  end
  if ~isempty(tfaces)
    TF = vertcat(tfaces{:});
  else
    TF = zeros(0, 3);
  end
  if ~isempty(norms)
    N = vertcat(norms{:});
  else
    N = zeros(0, 3);
  end
  if ~isempty(nfaces)
    NF = vertcat(nfaces{:});
  else
    NF = zeros(0, 3);
  end
  SI = zeros(0, 2);
end
