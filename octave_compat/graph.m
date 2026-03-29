function G = graph(s, t, w)
% GRAPH  Octave compatibility shim for MATLAB's graph object.
%   G = graph(s, t, w) creates a weighted undirected graph.
  if nargin < 3
    w = ones(size(s));
  end
  if isempty(s) && isempty(t)
    n = 0;
  else
    n = round(double(max([s(:); t(:)])));
  end
  G.s = round(double(s(:)));
  G.t = round(double(t(:)));
  G.w = double(w(:));
  G.n = n;
  % Build adjacency matrix (symmetric, using minimum weight for duplicates)
  G.adj = sparse(G.s, G.t, G.w, n, n) + sparse(G.t, G.s, G.w, n, n);
end
