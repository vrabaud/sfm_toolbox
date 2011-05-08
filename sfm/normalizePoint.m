function [ xHat T ] = normalizePoint( x, method )
% Normalize 2D/3D points for SFM
%
% The points can contain NaN and the entries are ignored in the case of
% isotropic scaling.
%
% USAGE
%  [ xHat T ] = normalizePoint(x,method)
%
% INPUTS
%  x       - [ nCoord x nPoint x nSet ] point coordinates (nSet can be 1)
%  method  - n (divide by the nth coordinate)
%            If Inf/-Inf, isotropic scaling (mean=sqrt(2) and centered
%            on 0).
%            If negative, the coordinates are sent back homogeneous.
%
% OUTPUTS
%  xHat    - the normalized points [ nCoord x nPoint x nSet ] or +1
%  T       - if method==Inf, it returns the transformation such that
%            xHat=T*x  (T is [ nCoord x nCoord x nSet ])
%
% EXAMPLE
%
% See also
%
% Vincent's Structure From Motion Toolbox      Version 3.1
% Copyright (C) 2008-2011 Vincent Rabaud.  [vrabaud-at-cs.ucsd.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the GPL [see external/gpl.txt]

xHat=x; nCoord = size(x,1); nSet = size(x,3);

if nargin<2 || isempty(method); method = nCoord; end

xIsnan=isnan(x); 
if any(xIsnan(:)); xIsnan=sum(xIsnan,1)>0; else xIsnan=[]; end

if abs(method)==Inf
  if ~isempty(xIsnan)
    c=zeros(nCoord, 1, nSet);
    scale=zeros(1, 1, nSet);
    for i=1:nSet
      % isotropic scaling
      % Get the centroid
      x_mask = x(:,~xIsnan(1,:,i),i);
      c(:,1,i) = mean( x_mask, 2 );
      % Get the scale factor
      scale(1,1,i) = sqrt(2)/mean(sqrt(sum( bsxfun(@minus, x_mask, ...
         c(:,1,i)).^2, 1)),2);
    end
  else
    % isotropic scaling
    % Get the centroid
    c = mean( x, 2 );
    % Get the scale factor
    scale = sqrt(2)./mean(sqrt(sum( bsxfun(@minus, x, c).^2, 1)),2);
  end


  T = repmat(eye(nCoord),[1,1,nSet]);
  T(:,end+1,:) = -c; T = bsxfun(@times,scale,T);
  
  xHat = bsxfun(@times, scale, bsxfun(@minus, x, c));
  
  if method<0; xHat(end+1,:,:) = 1; T(end+1,end,:)=1; end
else
  if abs(method)>nCoord
    if method<0; xHat( nCoord + 1, :, : ) = 1; end
  else
    finite = abs(x( abs(method), :))>eps;
    notMethod = 1 : nCoord; notMethod( abs(method) ) = [];
    xHat( notMethod, finite ) = bsxfun(@rdivide, x( notMethod, finite ),...
      x( abs(method), finite ));
    if method==nCoord
      xHat( method, :, : ) = [];
    else
      xHat( abs(method), : ) = finite;
    end
  end
end

if ~isempty(xIsnan)
  for i=1:nSet
    xHat(:,xIsnan(1,:,i),i)=NaN;
  end
end
