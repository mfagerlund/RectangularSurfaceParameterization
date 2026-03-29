function [bins, binsizes] = conncomp(G)
% CONNCOMP  Octave compatibility shim for MATLAB's conncomp.
%   [bins, binsizes] = conncomp(G) finds connected components.
  n = G.n;
  visited = false(n, 1);
  bins = zeros(n, 1);
  comp = 0;

  for start = 1:n
    if visited(start)
      continue;
    end
    comp = comp + 1;
    queue = start;
    visited(start) = true;
    bins(start) = comp;
    while ~isempty(queue)
      node = queue(1);
      queue(1) = [];
      neighbors = find(G.adj(node, :));
      for k = 1:length(neighbors)
        nb = neighbors(k);
        if ~visited(nb)
          visited(nb) = true;
          bins(nb) = comp;
          queue(end+1) = nb;
        end
      end
    end
  end

  binsizes = zeros(1, comp);
  for i = 1:comp
    binsizes(i) = sum(bins == i);
  end
end
