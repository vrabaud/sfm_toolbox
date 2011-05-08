function [grt, grR, grS] = msfmGradientSRt( anim, lamt, lamR, lamS )
% Compute the smoothness gradient in MSFM from Rabaud CVPR08
%
% Just check Rabaud's CVPR08 paper
% the notations follow Rabaud's PhD thesis: Appendix D3
%
% This function only deals with the part of the gradient
% with S, R and t, not with the dimensionality of S
%
% USAGE
%  [grt, grR, grS] = msfmGradientSRt( anim, lamt, lamR, lamS )
%
% INPUTS
%  anim          - Animation object (help Animation for details)
%  lamt          - weight of the translation smoothness error
%  lamR          - weight of the rotation smoothness error
%  lamS          - weight of the 3D shape smoothness error
%
% OUTPUTS
%  grt           - gradient in translation
%  grR           - gradient in rotation ([ 4 x nFrame], quaternionwise
%  grS           - gradient in shape
%
% EXAMPLE
%
% See also MSFMGRADIENTLSML
%
% Vincent's Structure From Motion Toolbox      Version 3.0
% Copyright (C) 2008-2011 Vincent Rabaud.  [vrabaud-at-cs.ucsd.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the GPL [see external/gpl.txt]

nFrame = anim.nFrame;

% compute ti*1n
ti1n = reshape(anim.t(1:2,:), 2, 1, nFrame);

% compute Ri*Si+ti*1n-Wi
firstParenthesis = 2*bsxfun(@plus, multiTimes( anim.R(1:2,:,:), ...
  anim.S, 2 ) - anim.W, ti1n);

% gradient with respect to ti
grt = reshape(sum(firstParenthesis,2),2,nFrame) + ...
  (2*lamt)*diffPrevNext(anim.t(1:2,:,:));
grt(3,:) = 0;

% gradient with respect to Ri
dgUdU = multiTimes( firstParenthesis, anim.S, 2.2 ) + ...
  (2*lamR)*diffPrevNext(anim.R(1:2,:,:));

% Now, chain rule: dg(U)/da,b,c,d=tr((dg(U)/dU)'*dU/da,b,c,d)
% where (a,b,c,d) is the quaternion
grR = zeros(4,nFrame);
dUdabcd = msfmRotationDerivative(quaternion(anim.R));
for i = 1 : 4
  grR(i,:) = multiTimes( dgUdU, squeeze(dUdabcd(:,:,i,:)), 4.1 );
end

% gradient with respect to Si
grS = multiTimes( anim.R(1:2,:,:), firstParenthesis, 2.1 ) + ...
  (2*lamS)*diffPrevNext(anim.S);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function res = diffPrevNext(X)
% the first element is X2-X1, and then 2Xi-X_{i-1}-X_{i+1}
if length(size(X))==2
  res = [ X(:,1)-X(:,2), -diff(X,2,2), X(:,end)-X(:,end-1) ];
else
  res = zeros(size(X));
  res(:,:,1) = X(:,:,1)-X(:,:,2);
  res(:,:,2:end-1) = -diff(X,2,3);
  res(:,:,end) = X(:,:,end)-X(:,:,end-1);
end
end
