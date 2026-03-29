function [Tree, pred] = minspantree(G, varargin)
% MINSPANTREE  Octave compatibility shim for MATLAB's minspantree.
%   [Tree, pred] = minspantree(G, 'Type','tree', 'Root',root)
%   Uses Prim's algorithm.

  p = inputParser;
  addParameter(p, 'Type', 'tree');
  addParameter(p, 'Root', 1);
  addParameter(p, 'Method', 'dense');
  parse(p, varargin{:});
  root = p.Results.Root;

  n = G.n;
  inMST = false(n, 1);
  key = inf(n, 1);
  pred = zeros(1, n);
  key(root) = 0;

  for iter = 1:n
    % Pick minimum key vertex not in MST
    k = key;
    k(inMST) = inf;
    [minval, u] = min(k);
    if minval == inf
      % Disconnected component - for forest mode, find next unvisited
      unvisited = find(~inMST, 1);
      if isempty(unvisited)
        break;
      end
      u = unvisited;
      key(u) = 0;
    end
    inMST(u) = true;

    % Update neighbors
    neighbors = find(G.adj(u, :));
    for k_idx = 1:length(neighbors)
      v = neighbors(k_idx);
      w = G.adj(u, v);
      if ~inMST(v) && w < key(v)
        key(v) = w;
        pred(v) = u;
      end
    end
  end

  % Build tree graph from pred
  s = [];
  t = [];
  w = [];
  for i = 1:n
    if pred(i) > 0
      s(end+1) = pred(i);
      t(end+1) = i;
      w(end+1) = G.adj(pred(i), i);
    end
  end
  if isempty(s)
    Tree = graph([], [], []);
    Tree.n = n;
  else
    Tree = graph(s(:), t(:), w(:));
    Tree.n = n;
  end
end
