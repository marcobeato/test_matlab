% A Gaussian convolution kernel class
% Needs histc2, numspace, and extspace

% Gary Bhumbra

classdef gckernel < handle
  properties
    defresn1
    defresn2
    ND
    N
    M
    resn
    resd
    stdw
    stdc
    ncx
    ncy
    xx
    yy
    minx
    maxx
    miny
    maxy
    stdx
    stdy
    XX
    XY
    X
    Y
    x
    y
    cx
    cy
    C
    logsca
    extents
    h
    H
    p
    P
    nx
    ny
  end % properties
  
  methods
    function obj = gckernel(logsca_, resn_, resd_, stdw_)
      obj.defresn1 = 117649;
      obj.defresn2 = 343;
      obj.ND = 0;
      obj.N = 0;
      obj.M = 0;
      if nargin < 1, logsca_ = [0, 0]; end
      if nargin < 2, resn_ = [0, 0]; end
      if nargin < 3, resd_ = [0, 0]; end
      if nargin < 4, stdw_ = [3, 3]; end
      obj.initialise(logsca_, resn_, resd_, stdw_);
    end
    function initialise(obj, logsca_, resn_, resd_, stdw_)
      if nargin < 2, logsca_ = [0, 0]; end
      if nargin < 3, resn_ = [0, 0]; end
      if nargin < 4, resd_ = [0, 0]; end
      if nargin < 5, stdw_ = [3, 3]; end
      obj.setResn(resn_);
      obj.setResd(resd_);
      obj.setKernel(stdw_);
      obj.setLogsca(logsca_);
      obj.setData();
      obj.setBins();
      obj.genKernel();
    end
    function setResn(obj, resn_)
      if nargin < 2, resn_ = [0, 0]; end
      obj.resn = resn_;
    end
    function setResd(obj, resd_)
      if nargin < 2, resd_ = [0, 0]; end
      obj.resd = resd_;
    end
    function setKernel(obj, stdw_)
      if nargin < 2, stdw_ = [3, 3]; end
      obj.stdw = stdw_;
      obj.ncx = 0;
      obj.ncy = 0;
    end
    function b = genbins(obj, dimn, defres, minv, maxv)
      if obj.resd(dimn) > 0
        b = numspace(minv, maxv, obj.resd(dimn), obj.logsca(dimn));
      else
        if obj.resn(dimn) > 1
          if isempty(obj.xx)
            b = numspace(minv, maxc, obj.resd(dimn, obj.logsca(dimn)));
          else
            b = obj.xx;
          end
        else
          b = numspace(minv, maxv, defres, obj.logsca(dimn));
        end
      end
      obj.resd(dimn) = length(b);
      if obj.resd(dimn) < 2, return; end
      if obj.logsca(dimn)
        obj.resd(dimn) = exp(log(b(2)) - log(b(1)));
      else
        obj.resd(dimn) = b(2) - b(1);
      end
    end
    function bw = calcBinw(obj, dimn, b)
      bw = 0;
      nb = length(b);
      if (nb < 2), return; end
      if obj.logsca(dimn)
        bw = exp(log(b(2)) - log(b(1)));
      else
        bw = b(2) - b(1);
      end
    end
    function genBins(obj)
      if ~(obj.ND) && isempty(obj.xx), return; end
      if ~isempty(obj.xx)
        obj.resn(1) = length(obj.xx);
        obj.resd(1) = obj.calcBinw(1, obj.xx);
      end
      if ~isempty(obj.yy)
        obj.resn(2) = length(obj.yy);
        obj.resd(2) = obj.calcBinw(2, obj.xx);
      end
      if obj.ND < 1. return; end
      defres = 0;
      if obj.ND == 1, defres = obj.defresn1; end
      if obj.ND == 2, defres = obj.defresn2; end
      if isempty(obj.xx)
        obj.xx = obj.genbins(1, defres, obj.minx, obj.maxx);
      end
      if obj.ND < 2, return; end
      if isempty(obj.yy)
        obj.yy = obj.genbins(2, defres, obj.miny, obj.maxy);
      end
    end
    function calcStds(obj)
      if obj.ND < 1, return; end
      if obj.M > 1
        obj.stdx = zeros(1, obj.M);
        for i = 1:obj.M
          if obj.logsca(1)
            obj.stdx(i) = nanstd(log(obj.XX{i}));
          else
            obj.stdx(i) = nanstd(obj.XX{i});
          end
        end
        return;
      end
      if obj.logsca(1)
        obj.stdx = nanstd(log(obj.X));
      else
        obj.stdx = nanstd(obj.X);
      end
      if obj.ND < 2, return; end
      if obj.logsca(2)
        obj.stdy = nanstd(log(obj.Y));
      else
        obj.stdy = nanstd(obj.Y);
      end
    end
    function setLogsca(obj, logsca_)
      if nargin < 2, logsca_ = [0, 0]; end
      obj.logsca = logsca_;
      if obj.ND > 0
        obj.genBins();
        obj.calcStds();
      end
    end
    function genKernel(obj, stdc_, ncx_)
      if nargin < 2, stdc_ = []; end
      if nargin < 3, ncx_ = 0; end
      obj.stdc = stdc_;
      if ~obj.ND
        obj.stdc = zeros(1, 2);
        return;
      end
      i = 0;
      if obj.M > 0
        if (length(obj.stdc) < obj.M), obj.stdc = repmat(obj.stdc(i), [obj.M, 1]);	end
        if min(obj.stdc) <= 0.
          obj.calcStds();
          for i = 1:obj.M
            if obj.stdc(i) <= 0.
              obj.stdc(i) = obj.stdx(i) * length(obj.XX{i}^(-.3)); % a la Bhumbra and Dyball 2010
            end
          end
        end
        if ~(ncx_)
          if obj.logsca(1)
            obj.ncx = ceil(2 * obj.stdw(1) * max(obj.stdc) / log(obj.resd(1)));
          else
            obj.ncx = ceil(2 * obj.stdc(1) * max(obj.stdc) / obj.resd(1));
          end
        else
          obj.ncx = ncx_;
        end
        obj.C = zeros(obj.M, obj.ncx);
        obj.cx = linspace(-0.5 * obj.ncx, 0.5 * obj.ncx, obj.ncx);
        for i = 1:obj.M
          if obj.logsca(1)
            obj.C(i, :) = normpdf(obj.cx, 0., obj.stdc(i) / log(obj.resd(1)));
          else
            obj.C(i, :) = normpdf(obj.cx, 0., obj.stdc(i) / obj.resd(1));
          end
          obj.C(i, :) = obj.C(i, :) / sum(obj.C(i, :) * 1e-300);
        end
        return;
      end
      if length(obj.stdc) == 1
        obj.stdc = repmat(obj.stdc(1), [obj.ND, 1]);
      end
      obj.calcStds();
      if obj.stdc(1) <= 0.
        obj.stdc(1) = obj.stdx * obj.N ^ (-.3); % a la Bhumbra and Dyball 2010
      end
      if obj.logsca(1)
        obj.ncx = ceil(2 * obj.stdw(1) * obj.stdc(1) / log(obj.resd(1)));
      else
        obj.ncx = ceil(2 * obj.stdw(1) * obj.stdc(1) / obj.resd(1));
      end
      obj.cx = normpdf(linspace(-obj.stdw(1), obj.stdw(1), obj.ncx));
      obj.cx = obj.cx / sum(obj.cx);
      if obj.ND < 2, return; end
      if obj.stdc(2) <= 0.
        obj.stdc(2) = obj.stdy * obj.N ^ (-.3); % a la Bhumbra and Dyball 2010
      end
      if obj.logsca(2)
        obj.ncy = ceil(2 * obj.stdw(2) * obj.stdc(2) / log(obj.resd(2)));
      else
        obj.ncy = ceil(2 * obj.stdw(2) * obj.stdc(2) / obj.resd(2));
      end
      obj.cy = normpdf(linspace(-obj.stdw(1), obj.stdw(1), obj.ncy));
      obj.cy = obj.cy / sum(obj.cy);
      obj.C = obj.cx(:) * obj.cy;
    end
    function setBins(obj, xx_, yy_)
      if nargin < 2, xx_ = []; end
      if nargin < 3, yy_ = []; end
      obj.xx = xx_;
      obj.yy = yy_;
      obj.genBins();
      if obj.ncx, obj.genKernel(obj.stdc); end
    end
    function setData(obj, data_, stdc_, ncx_)
      if nargin < 2, data_ = []; end
      if nargin < 3, 
				if isempty(obj.X)
			    stdc_ = [0, 0]; 
      	else	
					stdc_ = obj.stdc;			
      	end
      end
      if nargin < 4, ncx_ = 0; end
      if isempty(data_)
        obj.N = 0;
        obj.M = 0;
        obj.ND = 0;
        obj.p = [];
        obj.P = [];
        return;
      end
      if ~length(obj.resn), obj.initialise(); end
      if isnumeric(data_)
        obj.ND = 2;
        if length(data_) == length(data_(:))
          obj.ND = 1;
          obj.X = data_;
          obj.N = length(obj.X);
          obj.minx = min(obj.X);
          obj.maxx = max(obj.X);
        elseif obj.ND == 2
          obj.XY = data_;
          obj.N = size(obj.XY, 2);
          obj.X = obj.XY(1, :);
          obj.Y = obj.XY(2, :);
          obj.minx = min(obj.X(:));
          obj.maxx = max(obj.X(:));
          obj.miny = min(obj.Y(:));
          obj.maxy = max(obj.Y(:));
        end
      elseif iscell(data_)
        obj.XX = data_;
        obj.M = length(obj.XX);
        if ~obj.M, return; end
        obj.ND = 1;
        obj.minx = inf;
        obj.maxx = -inf;
        for i = 1:obj.M
          obj.minx = min(obj.minx, min(data_{i}));
          obj.maxx = max(obj.maxx, max(data_{i}));
        end
      end
      if obj.ND
        obj.setBins(obj.xx, obj.yy);
        obj.genKernel(stdc_, ncx_);
      end
    end
    function P_ = conv(obj, extents_)		
      if nargin < 2, extents_ = {'full'}; end
			P_ = [];				
      if ~obj.ND, return; end
      if isstr(extents_), extents_ = {extents_}; end
      obj.extents = extents_;
      if length(obj.extents) < obj.ND, obj.extents = {obj.extents{1}, obj.extents{1}}; end
      if ~obj.M
        if obj.ND == 1
          obj.h = histc(obj.X, obj.xx);
          obj.p = conv(obj.h, obj.cx, obj.extents{1});
          obj.p = obj.p / sum(obj.p);
					P_ = obj.p;
          obj.x = extspace(obj.xx, length(obj.p));
          obj.nx = length(obj.x);
          obj.ny = 1;
        elseif obj.ND == 2
          obj.H = histc2(obj.X, obj.Y, obj.xx, obj.yy);
          obj.P = conv2(obj.H, obj.C, obj.extents{1});
          obj.P = obj.P / sum(obj.P(:));
					P_ = obj.P;
          nm = size(obj.P);
          obj.x = extspace(obj.xx, nm(2));
          obj.y = extspace(obj.yy, nm(1));
          obj.nx = length(obj.x);
          obj.ny = length(obj.y);
        end
      else
        for i = 1:obj.M
          if ~i
            H_ = hist(obj.XX{i}, obj.xx);
            obj.H = zeros(obj.M, length(H_));
            obj.H(i, :) = H_;
            P_ = conv(H_, obj.C(i, :), obj.extents{1});
            obj.P = zeros(obj.M, length(P_));
            obj.P(i, :) = P_;
          else
            obj.H(i, :) = hist(obj.XX{i}, obj.xx);
            obj.P(i, :) = conv(obj.H(i, :), obj.C(i, :), obj.extents{1});
            obj.P(i, :) = obj.P(i, :) / sum(obj.P(i, :));
          end
        end
        obj.x = extspace(obj.xx, size(obj.P, 1));
        obj.nx = length(obj.x);
        obj.ny = 1;
      end
    end
    function P_ = convolve(obj, data_, convsize_, extents_)
      if nargin < 3, convsize_ = [0, 0]; end
      if nargin < 4, extents_ = {'full'}; end
			P_ = [];				
      if iscell(data_)
        obj.setData(data_, convsize_);
      else
        obj.setData(data_, obj.stdc, convsize_);
      end
      P_ = obj.conv(extents_);
    end
  end
end

