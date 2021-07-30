function [fxy, xn, yn] = hist2c(xval, yval, xbin, ybin);
% function [fxy, xn, yn] = hist2c(xval, yval, xbin, ybin);
%
% Generates a two-dimensional histogram matrix of paired
% values xval and yval divided by rectangular bins whose
% edges are given by xbin and ybin. The histogram heights
% are given by fxy whereas the bin references for each of
% the values are given by xn and yn. If ybin is omitted,
% xbin is used for both xval and yval.

% Revised on 20th-vi-MMII by Gary Bhumbra.

if nargin<4
    ybin = xbin;
end

xval = xval(:);
yval = yval(:);

nxval = length(xval);
nyval = length(yval);

if nxval ~= nyval;
    error('X & Y arrays are not of the same length');
    % nval = min(nxval,nyval);
    % xval = xval(1:nval);
    % yval = yval(1:nval);
else
    nval = nxval;
end

xbin = xbin(:);
ybin = ybin(:);

nxbin = length(xbin);
nybin = length(ybin);

[xh, xn] = histc(xval, xbin);
[yh, yn] = histc(yval, ybin);

nnr = repmat(1:nxbin, [nval 1]);
nnl = (0.5:1:nybin+0.5)';

xnr = repmat(xn, [1 nxbin]);
ynr = repmat(yn, [1 nxbin]);

xy = ynr.*logical(xnr==nnr);
fxy = zeros(nybin+1, nxbin);

for f = 1:nxbin;
    fxy(:, f) = histc(xy(:, f), nnl);
end

fxy = fxy(1:end-1, :);
