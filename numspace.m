function x = numspace(minx, maxx, n, uselog)
% function x = numspace(minx, maxx, n, uselog)

% Gary Bhumbra

if nargin < 1, minx = 0; end
if nargin < 2, maxx = 1; end
if nargin < 3, n = 100; end
if nargin < 4, uselog = false; end

if uselog
	x = logspace(log10(minx), log10(maxx), n);
else
  x = linspace(minx, maxx, n);
end

