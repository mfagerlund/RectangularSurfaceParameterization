function x = quadprog(H, f, Ain, bin, Aeq, beq, lb, ub, x0)
% QUADPROG  Wrapper that handles rank-deficient equality constraints.
%   Octave's quadprog requires full row rank Aeq. This removes redundant rows.

  if nargin < 9, x0 = []; end
  if nargin < 8, ub = []; end
  if nargin < 7, lb = []; end
  if nargin < 6, beq = []; end
  if nargin < 5, Aeq = []; end
  if nargin < 4, bin = []; end
  if nargin < 3, Ain = []; end

  if ~isempty(Aeq) && size(Aeq, 1) > 1
    % Remove near-zero rows
    row_norms = sqrt(sum(Aeq.^2, 2));
    keep = row_norms > 1e-14;
    Aeq = Aeq(keep, :);
    beq = beq(keep);

    % Use SVD to find rank and select independent rows
    if size(Aeq, 1) > 1
      [U, S, V] = svd(full(Aeq), 'econ');
      s = diag(S);
      tol = max(size(Aeq)) * eps(s(1));
      r = sum(s > tol);
      if r < size(Aeq, 1)
        % Project to independent set: use first r rows of U'*Aeq
        Aeq_new = U(:,1:r)' * Aeq;
        beq_new = U(:,1:r)' * beq;
        Aeq = Aeq_new;
        beq = beq_new;
      end
    end
  end

  % Call the real quadprog from optim package via __octave_config_info__
  % We need to avoid recursion - call the builtin directly
  real_qp = @__qp__;

  % Convert to __qp__ format: min 0.5*x'*H*x + f'*x
  % subject to: Aeq*x = beq, Ain*x <= bin, lb <= x <= ub
  n = length(f);
  if isempty(lb), lb = -inf(n,1); end
  if isempty(ub), ub = inf(n,1); end
  if isempty(x0), x0 = zeros(n,1); end
  if isempty(Ain), Ain = zeros(0,n); bin = zeros(0,1); end
  if isempty(Aeq), Aeq = zeros(0,n); beq = zeros(0,1); end

  [x, ~, info] = qp(x0, H, f, Aeq, beq, lb, ub, [], Ain, bin);
  if info.info ~= 0
    warning('quadprog: solver returned info=%d', info.info);
  end
end
