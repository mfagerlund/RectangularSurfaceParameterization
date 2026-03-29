function P = shortestpath(G, source, target)
% SHORTESTPATH  Octave compatibility shim for MATLAB's shortestpath.
%   P = shortestpath(G, source, target) finds shortest path using Dijkstra.
  n = G.n;
  dist = inf(n, 1);
  prev = zeros(n, 1);
  dist(source) = 0;
  visited = false(n, 1);

  for iter = 1:n
    % Find unvisited node with smallest distance
    d = dist;
    d(visited) = inf;
    [~, u] = min(d);
    if dist(u) == inf
      break;
    end
    if u == target
      break;
    end
    visited(u) = true;

    % Relax neighbors
    neighbors = find(G.adj(u, :));
    for k = 1:length(neighbors)
      v = neighbors(k);
      if ~visited(v)
        alt = dist(u) + G.adj(u, v);
        if alt < dist(v)
          dist(v) = alt;
          prev(v) = u;
        end
      end
    end
  end

  % Reconstruct path
  P = target;
  u = target;
  while prev(u) ~= 0
    u = prev(u);
    P = [u, P];
  end
end
