function out = quaternion( in, doSmooth )
% Converts between rotation matrix and quaternion representation
%
% The formulas were taken from wikipedia:
% http://en.wikipedia.org/wiki/Rotation_matrix
%
% USAGE
%  Q = quaternion( R )
%  R = quaternion( Q )
%
% INPUTS
%  in       - either 3x3xN set of rotation matrices or 4xN or Nx4
%             set of quaternions
%  doSmooth - if set to true, the output quaternions are opposed if
%             necessary to get smoother data
%
% OUTPUTS
%  out      - 4xN set of quaternion or 3x3xN set of rotation matrices
%
% EXAMPLE
%  R = zeros(3,3,1000);
%  for i=1:1000
%    R(:,:,i) = rotationMatrix( rand(1, 3) );
%  end
%  norm( sum((quaternion(quaternion(R)) - R).^2,3), 'fro' )
%
% See also
%
% Vincent's Structure From Motion Toolbox      Version 3.0
% Copyright (C) 2008-2011 Vincent Rabaud.  [vrabaud-at-cs.ucsd.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the GPL [see external/gpl.txt]

if nargin<=1; doSmooth = false; end

if size(in,1)==3
  nFrame = size(in,3);
else
  if size(in,1)~=4; in = in'; end
  nFrame = size( in, 2 );
end

switch size(in,1)
  case 3
    out = zeros( 4, nFrame );
    
    r = zeros(1,nFrame); s=r; w=r; x=r; y=r; z=r;
    
    Qxx = squeeze(in(1,1,:))'; Qyx = squeeze(in(2,1,:))'; Qzx = squeeze(in(3,1,:))';
    Qxy = squeeze(in(1,2,:))'; Qyy = squeeze(in(2,2,:))'; Qzy = squeeze(in(3,2,:))';
    Qxz = squeeze(in(1,3,:))'; Qyz = squeeze(in(2,3,:))'; Qzz = squeeze(in(3,3,:))';
    
    t = Qxx+Qyy+Qzz+1.0;
    
    ind=find(t>=1e-3);
    r(ind) = sqrt(t(ind));
    s(ind) = 0.5./r(ind);
    w(ind) = 0.5.*r(ind);
    x(ind) = (Qzy(ind)-Qyz(ind)).*s(ind);
    y(ind) = (Qxz(ind)-Qzx(ind)).*s(ind);
    z(ind) = (Qyx(ind)-Qxy(ind)).*s(ind);
    
    ind=find( t<1e-3 & Qxx >= Qyy & Qxx >= Qzz);
    r(ind) = sqrt( 1.0 + Qxx(ind) - Qyy(ind) - Qzz(ind)); s(ind) = 0.5./r(ind);
    w(ind) = (Qzy(ind) - Qyz(ind) ).*s(ind); x(ind) = 0.5 * r(ind);
    y(ind) = (Qxy(ind) + Qyx(ind) ).*s(ind); z(ind) = (Qxz(ind) + Qzx(ind) ).*s(ind);
    
    ind=find( t<1e-3 & Qyy >= Qxx & Qyy >= Qzz);
    r(ind) = sqrt( 1.0 + Qyy(ind) - Qxx(ind) - Qzz(ind)); s(ind) = 0.5./r(ind);
    w(ind) = (Qxz(ind) - Qzx(ind) ).*s(ind); x(ind) = (Qxy(ind) + Qyx(ind) ).*s(ind);
    y(ind) = 0.5 * r(ind); z(ind) = (Qyz(ind) + Qzy(ind) ).*s(ind);
    
    ind=find( t<1e-3 & Qzz >= Qxx & Qzz >= Qyy );
    r(ind) = sqrt( 1.0 + Qzz(ind) - Qxx(ind) - Qyy(ind) ); s(ind) = 0.5./r(ind);
    w(ind) = (Qyx(ind) - Qxy(ind) ).*s(ind); x(ind) = (Qxz(ind) + Qzx(ind) ).*s(ind);
    y(ind) = (Qyz(ind) + Qzy(ind) ).*s(ind); z(ind) = 0.5 .* r(ind);
    
    out = [ w; x; y; z ];
    
    % Normalize
    out = bsxfun(@rdivide, out, sqrt(sum(out.^2,1)));
    
    % Smooth if requested
    if doSmooth
      for i=2:nFrame
        if norm(out(:,i)-out(:,i-1)) > norm(-out(:,i)-out(:,i-1))
          out(:,i) = -out(:,i);
        end
      end
    end
  otherwise
    s = zeros(1,nFrame);
    w = in(1,:); x = in(2,:); y = in(3,:); z = in(4,:);
    
    Nq = w.*w + x.*x + y.*y + z.*z;
    s(Nq>0)=2./Nq(Nq>0);
    
    X = x.*s;
    Y = y.*s;
    Z = z.*s;
    wX = w.*X;
    wY = w.*Y;
    wZ = w.*Z;
    xX = x.*X;
    xY = x.*Y;
    xZ = x.*Z;
    yY = y.*Y;
    yZ = y.*Z;
    zZ = z.*Z;
    
    out = reshape( [ 1.0-(yY+zZ); xY+wZ; xZ-wY; xY-wZ; 1.0-(xX+zZ); yZ+wX; xZ+wY; yZ-wX; 1.0-(xX+yY) ], 3, 3, nFrame );
end
