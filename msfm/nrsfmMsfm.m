function anim = nrsfmMsfm( W, d, animGT )
% Compute orthographic non-rigid structure from motion using MSFM from Rabaud
%
% Just check Rabaud's CVPR08 paper
%
% USAGE
%  anim = nrsfmMSFM( W )
%
% INPUTS
%  W              - [ 2 x nPoint x nFrame ] set of 2D points
%  d              - the dimensionality of the embedding
%  animGT         - the ground truth Animation. This is obviously not
%                   necessary but it will display how the 3D error changes
%                   iteration after iteration. I will remove that once the
%                   code is faster.
%
% OUTPUTS
%  anim           - Animation object (help Animation for details)
%
% EXAMPLE
%
% See also
%
% Vincent's Structure From Motion Toolbox      Version 3.0
% Copyright (C) 2008-2011 Vincent Rabaud.  [vrabaud-at-cs.ucsd.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the GPL [see external/gpl.txt]

doAnimGT=(nargin>=3);

nFrame = size(W,3); nPoint = size(W,2);
lamt=1*0; lamR=1/100*0; lamS=1*0;
anim = Animation(); anim.W = W;

% Compute pairwise similarities
Sim2=computeAnimSimilarity(anim,2);

% Compute triadic similarities
%Sim3=computeAnimSimilarity(anim,3,Sim2);

% Compute the rigid shape chain
%[animRSC,basis,ind] = msfmComputeRsc( anim, Sim2, Sim3);
[animRSC,basis,ind] = msfmComputeRsc( W, Sim2 );

% perform one iteration of gradient descent to spread out the clusters
% but by only using the error criterion, not the dimensionality one,
% as neighborhoods are just points right now, they are not flat enough
error = errorMsfm(animRSC, lamt, lamR, lamS);

if doAnimGT
  disp(error)
  disp(animRSC.computeError( 'animGT', animGT))
  fflush(stdout);
end
[grt, grR, grS] = msfmGradientSRt( animRSC, lamt, lamR, lamS );
anim = oneGradientDescentIter(animRSC, zeros(3,nPoint,nFrame), grt, ...
  grR, grS, lamt, lamR, lamS);
error = errorMsfm(anim, lamt, lamR, lamS);
if doAnimGT
  disp(error)
  disp(anim.computeError( 'animGT', animGT))
  fflush(stdout);
end
% perform gradient descent to reach the local optimum
disp('starting final optimization');
err=error(5); errPrev = 2*err;
while (errPrev-err)>0.001*errPrev
  errPrev=err;
  [ grt, grR, grS ] = msfmGradientSRt( anim, lamt, lamR, lamS );
  [ grSOrtho, H ] = msfmGradientLsml( anim, d );
  
  % only keep the component of grS that is on the manifold
  % as sticking to a low-dimensional manifold is the priority
  grSTangent = grS - reshape( bsxfun( @times, sum( bsxfun(@times, ...
    reshape(grS,3*nPoint,1,nFrame), H), 1 ), H ), 3, nPoint, nFrame );
  
  % perform line search
  anim = oneGradientDescentIter(anim, grSOrtho, grt, grR, ...
    grSTangent, lamt, lamR, lamS);
  err=errorMsfm(anim, lamt, lamR, lamS)
  err=err(5);
  if doAnimGT, disp(anim.computeError('animGT', animGT)); end
  fflush(stdout);
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function animBest = oneGradientDescentIter(anim, grSOrtho, grt, grR, ...
  grSTangent, lamt, lamR, lamS)
% cache some data
errorOld = errorMsfm(anim, lamt, lamR, lamS);
animOld = anim;
oldRQuaternion = quaternion(animOld.R);

% make grSTangent of the same norm as grSOrtho for each point
nPoint=size(grSTangent,2); nFrame=size(grSTangent,3);

% try with different length of grSTangent
animBest = anim; errBest = errorOld(5);
isDone = false; grSTangentOri = grSTangent;
for lam1 = 0 : 10
  for lam2 = 0 : 10
    grS = 1/2^lam2*grSTangent + grSOrtho;
    lam = 1/2^lam1;
    anim.t = animOld.t - lam*grt;
    anim.R = quaternion(oldRQuaternion - lam*grR);
    anim.S = animOld.S - lam*grS;
    error = errorMsfm(anim,lamt,lamR,lamS);
    if error(5) < errBest, animBest=anim; isDone = true; break; end
  end
end
end

function errTot = errorMsfm(anim, lamt, lamR, lamS)
errTot = zeros(1,5);

% reprojection error
err = anim.computeError();
errTot(1) = err(1);

% error in t
errTot(2) = lamt*norm( diff(anim.t,1,2), 'fro' )^2;

% error in R
errTot(3) = lamR*norm( reshape( diff(anim.R,1,3), 1, [] ))^2;

% error in S
errTot(4) = lamS*norm( reshape( diff(anim.S,1,3), 1, [] ))^2;

% total error
errTot(5) = sum( errTot );
end
