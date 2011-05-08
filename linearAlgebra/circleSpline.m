function [ p dp ] = circleSpline( nPoint, doDisplay, nIter )
% Returns a smooth random paths on a sphere
%
% This is an implementation of Fair, G - and C -Continuous Circle Splines for the Interpolation of Sparse Data JCAD 05, by Sequin and Lee
% Explained in: www.cs.berkeley.edu/~sequin/PAPERS/JCAD05_CircleSplines.pdf
% www.cs.berkeley.edu/~sequin/TALKS/SIGG03_Csplines.ppt
%
% USAGE
%  pu= circleSpline( nPoint )
%
% INPUTS
%  nPoint       - number of points on the unitsphere
%  doDisplay    - [false] plot the final data
%  nIter        - [4] number of times the path is repeated
%
% OUTPUTS
%  p         - [ 3 dimSSpace  x nFrame ] generated 3D points
%  p         - [ 3 dimSSpace  x nFrame ] corresponding derivatives
%
% EXAMPLE
%
% See also
%
% Vincent's Structure From Motion Toolbox      Version 2.11
% Copyright (C) 2008-2011 Vincent Rabaud.  [vrabaud-at-cs.ucsd.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the GPL [see external/gpl.txt]

if nargin<=1; doDisplay = false; end
if nargin<=2 || nIter<=1; nIter = 4; end

p = zeros(3,nIter*4);

while 1; ang = [ 0 rand(1,3) ]*2*pi; if any(abs(diff(ang))>pi/2); break; end; end
% make sure 2 heights are positive and 2 negative
while 1; height = 2*rand(1,4)-1; if sum(height>0)==2 && max(height)>0.75 && min(height)<-0.75; break; end; end
k = 1;
% Generate several key points (with a certain symmetry)
for i = 0 : 2*pi/nIter : 2*pi
  for j = 1 : 4
    p(:,k) = [cos(i+ang(j));sin(i+ang(j));height(j)]; k = k + 1;
  end
end
p(:,end) = []; % just to have it loop
p(1:2,:) = p(1:2,:)./repmat(sqrt(sum(p(1:2,:).^2,1))./sqrt(1-p(3,:).^2),[2 1]);
%  minmax(sum(p.^2,1))

%  p = [ 1 -0.5 0; 0 0 0; 3 0 0; 2 0.5 0 ]';
%  p = [ 0 0 -1; 0 0 0; 1 0 0; 1 -1 0 ]';
%  p = [ 0 0.5 0; 0.5 0 0; 1.5 0.5 0; 1.25 1.2 0 ]';
%  p = [ 0 0 0; 0 1 0; 1 0 0; 1 1 0 ]';
% Treat consecutive overlaping triplets and compute points in between
%  clf; hold on;
k = 1; pu = zeros(3,1000);
for i = 2 : size(p,2) - 2
  % Get the two tangent vectors
  a = p(:,i) - p(:,i-1); a = a/norm(a);
  b = p(:,i+1) - p(:,i); bNorm = norm(b); b = b/bNorm;
  c = p(:,i+1) - p(:,i-1); c = c/norm(c);
  d = p(:,i+2) - p(:,i+1); d = d/norm(d);
  e = p(:,i+2) - p(:,i); e = e/norm(e);
  
  tau = acos( [ a'*c e'*d ] ); % based on beta
  axist = [ cross(b,a) cross(b,d) ];
  
  % Compute the right tangent (and not like described in the paper)
  % as the taus are good up to an ambiguity
  t = zeros(3,2);
  for j = 1 : 2
    tTmp = [ rotationMatrix(axist(:,j),tau(j))*b, rotationMatrix(axist(:,j),-tau(j))*b ];
    % Check which one is more tangent
    if j==1
      alpha = acos( (-b)'*(-c) );
      alphaTmp = acos( a'*tTmp );
    else
      alpha = acos( b'*e );
      alphaTmp = acos( (-d)'*(-tTmp) );
    end
    [ disc ind ] = min( abs( alphaTmp - alpha ) );
    t(:,j) = tTmp(:,ind);
    if ind==2; tau(j) = -tau(j); end
  end
  
  % get the mirrored t(:,2)
  t2Orthob = t(:,2) - (t(:,2)'*b)*b;
  t(:,3) = t(:,2) - 2*t2Orthob;
  tau(2) = -tau(2);
  
  axisSwivel = cross(t(:,1),t(:,3));
  tauSwivelMax = acos(t(:,1)'*t(:,3));
  if pi-tauSwivelMax<0.01; continue; end
  % if p(:,i+1) is not on the right side
  if (t(:,1)+t(:,3))'*(p(:,i+1)-p(:,i))<0
    tauSwivelMax = -2*pi+tauSwivelMax;
  end
  
  if norm(rotationMatrix(axisSwivel,tauSwivelMax)*t(:,1)-t(:,3)) > norm(rotationMatrix(-axisSwivel,tauSwivelMax)*t(:,1)-t(:,3))
    axisSwivel = - axisSwivel;
  end
  for u = 0 : 0.01 : 1
    % Reflect tau(2)
    %      tauu = tau(1)*cos(u*pi/2)^2 + tau(2)*sin(u*pi/2)^2;
    tauSwivel = 0*cos(u*pi/2)^2 + tauSwivelMax*sin(u*pi/2)^2;
    tu = rotationMatrix(axisSwivel,tauSwivel)*t(:,1);
    tauu = acos(tu'*b);
    
    %      tmp1 = [ 0 1 0 ] + tu'; tmp1 = [ 0 1 0; tmp1 ];
    %      line( tmp1(:,1), tmp1(:,2), tmp1(:,3) );
    
    fu = bNorm*sin(u*tauu)/sin(tauu);
    if abs(tauu)<1e-3; fu = bNorm*u*sinc(u*tauu)/sinc(tauu); end
    
    phiu = (1-u)*tauu;
    %  [ fu phiu ]
    % Compute the P(u) point
    if norm(cross(tu,b))<1e-2; continue; end
    puTmp = fu*(rotationMatrix(cross(tu,b),-abs(phiu))*b);
    
    %      puTmp = fu*( cos(phiu)*b + sin(phiu)*tuOrtho );
    pu(:,k) = p(:,i) + puTmp;
    k = k + 1;
    if k>=size(pu,2); pu = [ pu zeros(3,1000) ]; end
  end
end
pu = pu(:,1:k-1);
dist = sqrt(sum(diff(pu,1,2).^2,1));
distCum = cumsum(dist);
len = sum(dist);

% Get the good subpoints
p = zeros( 3, nPoint ); dp = p;
kp = 2;
p(:,1) = pu(:,1);
dp(:,1) = pu(:,2)-pu(:,1);
for kpu = 2 : length(distCum)-1
  if distCum(kpu)<(kp-1)*len/(nPoint-1) && distCum(kpu+1)>(kp-1)*len/(nPoint-1)
    [ disc ind ] = min( abs(distCum(kpu:kpu+1)-(kp-1)*len/(nPoint-1)) );
    p(:,kp) = pu(:,kpu-1+ind); dp(:,kp) = pu(:,kpu-1+ind+1) - pu(:,kpu-1+ind-1);
    kp = kp + 1;
  end
end
p(:,end) = pu(:,end);
dp(:,end) = pu(:,end)-pu(:,end-1);

% Make sure the points are really on the unit sphere (remove small errors)
p = p./repmat(sqrt(sum(p.^2,1)),[3 1]);
if doDisplay
  clf; plot3(p(1,:),p(2,:),p(3,:),'.-'); hold on; axis equal;
  %    plot3(p(1,1),p(2,1),p(3,1),'sr');
  %    plot3(p(1,end),p(2,end),p(3,end),'*g','MarkerSize',20);
  %    plot3(p(1,:),p(2,:),p(3,:),'sr-');
end
if nargout>=2
  dp = dp./repmat(sqrt(sum(dp.^2,1)), [ 3 1 ] );
end
%    tmp1 = [ 0 1 0 ] + t(:,1)'; tmp1 = [ 0 1 0; tmp1 ];
%    tmp2 = [ 0 1 0 ] + t(:,3)'; tmp2 = [ 0 1 0; tmp2 ];
%  h1 = line( tmp1(:,1), tmp1(:,2), tmp1(:,3) );
%  h2 = line( tmp2(:,1), tmp2(:,2), tmp2(:,3) );
%  set(h1,'Color', 'g' ); set(h2,'Color', 'g' );
%  minmax(sum(p.^2,1))
