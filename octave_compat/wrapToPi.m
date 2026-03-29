function a = wrapToPi(a)
% WRAPTOPI  Wraps angles to [-pi, pi].
  a = mod(a + pi, 2*pi) - pi;
end
