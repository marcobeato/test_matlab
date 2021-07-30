function y = extspace(x, n, logsca)
if nargin < 2, n = 1; end
if nargin < 3, logsca = false; end
nx = length(x);				
minx = min(x);
maxx = max(x);
if logsca
  minx = log(minx);
  maxx = log(maxx);
end

midx = 0.5 * (minx + maxx);

if n == 1 || nx == 1
  y = repmat(midx, [1, max(n, nx)]); 
else
  widx = 0.5 * (maxx - minx);
	widy = widx * (n - 1) / (nx - 1);
  y = linspace(midx-widy, midx+widy, n);
end

if logsca
  y = exp(y);
end	
