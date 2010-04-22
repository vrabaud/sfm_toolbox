function anim = metricUpgrade(anim, isCalibrated)
% Perform an affine upgrade
%
% Given projection matrices and a 3D projective structure in anim,
% perform an affine upgrade
%
% USAGE
%   anim = affineUpgrade(anim, isCalibrated)
%
% INPUTS
%  anim       - Animation object with P and S filled
%
% OUTPUTS
%  anim      - Animation object such that 
%
% EXAMPLE
%
% See also
%
% Vincent's Structure From Motion Toolbox      Version 3.0
% Copyright (C) 2008-2010 Vincent Rabaud.  [vrabaud-at-cs.ucsd.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the GPL [see external/gpl.txt]

P=anim.P; S=anim.S; W=anim.W; nPoint=anim.nPoint; nFrame=anim.nFrame;

if ~anim.isProj
  if isCalibrated
    % Tomasi Kanade with the metric constraint
    G=zeros( 3*nFrame, 6 );
    for i=1:nFrame
      G(i,:) = g( P(1,:,i), P(1,:,i) );
      G(nFrame+i,:) = g( P(2,:,i), P(2,:,i) );
      G(2*nFrame+i,:) = g( P(1,:,i), P(2,:,i) );
    end
    % Solve for the square root matrix
    l = G\[ ones(2*nFrame,1); zeros(nFrame,1) ];
    L = [ l(1:3)'; l(2) l(4:5)'; l(3) l(5:6)' ];
    [ U S V ] = svd(L); S(S<0) = 0; Q = U*sqrt(S)*V';

    P(1:2,1:3,:)=multiTimes(P(1:2,1:3,:),Q,1);
    anim.R = zeros(3,3,nFrame);
    for i=1:nFrame
      anim.R(:,:,i) = rotationMatrix( rotationMatrix(P(1:2,1:3,i) ));
    end
    anim.t=reshape(P(1:2,4,:),2,nFrame); anim.t(3,:)=0;
    
    % get the optimal S
    anim.S=reshape(permute(P(1:2,1:3,:),[1,3,2]),2*nFrame,3)\...
      reshape(permute(bsxfun(@minus,W,P(1:2,4,:)),[1,3,2]),2*nFrame,nPoint);
  end
  return
end

% Chandraker IJCV 2009
% globally optimal metric upgrade

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function res=g(a,b)
% perform some simple product for Kanade factorization
res=[ a(1)*b(1) a(1)*b(2)+a(2)*b(1) a(1)*b(3)+a(3)*b(1) a(2)*b(2) ...
  a(2)*b(3)+a(3)*b(2) a(3)*b(3) ];
end
