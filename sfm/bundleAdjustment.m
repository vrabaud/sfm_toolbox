function [ anim info ] = bundleAdjustment( anim, varargin )
% Bundle Adjustment (using the excellent SBA)
%
% Can perform BA on full camera matrices (3 x 4) or K, R, t
%
% This can work for homographies too in case S0 is [ 2 x nPoint ]
%
% The calibration parameters are given in 3 or 5 by nFrame matrices.
% Each line corresponds to an element:
% Affine
%   K1  K2  0
%   0   K3  0
%   0   0   1
% Projective
%   K1  K2  K4
%   0   K3  K5
%   0   0   1
%
% USAGE
%  [ anim ] = bundleAdjustment( anim, sfmCase, varargin )
%
% INPUTS
%  anim      - Animation object
%  varargin   - list of paramaters in quotes alternating with their values
%       - sfmCase   - 1,2,3, Inf motion+struct, motion only, struct only,
%              NRSFM with shape basis
%       - 'KMask', [ 3 x 1 ] or [ 5 x 1 ] contains 1 when the calibration
%               parameter is fixed
%       - 'nItrSBA' number of BA iterations
%       - 'nFrameFixed' [1], number of frames (starting from the first)
%                       whose camera parameters are fixed
%
% OUTPUTS
%  anim      - optimized animation with optimized S, P, K, R, t
%
% EXAMPLE
%
% See also
%
% Vincent's Structure From Motion Toolbox      Version 3.0
% Copyright (C) 2009 Vincent Rabaud.  [vrabaud-at-cs.ucsd.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the GPL [see external/gpl.txt]

[ KMask nItrSBA covx nFrameFixed ] = ...
  getPrmDflt( varargin,{ 'KMask', [], 'nItrSBA', 500, ...
  'covx', [], 'nFrameFixed', 1 }, 1 );

% Figure out where the files are located
switch computer
  case 'PCWIN',
    projPath=[ '@' fileparts(mfilename('fullpath')) ...
      '/private/sba/sbaProjection.dll' ];
  case {'GLNX86','i686-pc-linux-gnu', 'GLNXA64', 'x86_64-pc-linux-gnu'}
    projPath=[ '@' fileparts(mfilename('fullpath')) ...
      '/private/sba/sbaProjection.so' ];
end

% set some parameters from the anim obejct
if isempty(anim.mask); WMask=ones(anim.nPoint, anim.nFrame);
else WMask=anim.mask; end

KAll=anim.K;

% Initialize the camera parameters with what already exists
% if there is nothing, use the full projection matrices
if ~isempty(anim.R) && ~isempty(anim.t)
  % create P0 for rotation and stuff
  if isempty(anim.K)
    K0=zeros(0,anim.nFrame); nK = 0; doK=0;
  else
    if isempty(KMask)
      if anim.isProj; KMask = ones(5,1); else KMask = ones(3,1); end
    end
    doK=1; K0=anim.K(~KMask,:); nK = length(find(~KMask));
  end
  
  % Deal with NRSFM
  switch size(anim.l,1)
    case anim.nBasis, %Xiao
      isFirstCoeff1=false;
    case anim.nBasis-1, %Torresani
      isFirstCoeff1=true;
    otherwise
      if ~isempty(anim.l); error('Problem with the dimension of l and SBasis');end
  end
  
  % Get the quaternion from the rotation matrices
  Q = quaternion(anim.R);
  
  %  sba(n, m, mcon, vmask, p0, cnp, pnp, x, covx, mnp, proj, projac, ...
  %itmax, verbose, opts, reftype, varargin)
  pnp=3;
  if anim.isProj % Projective camera
    P0 = [ Q; anim.t; K0 ]; cnp=7+nK;
  else % Affine camera
    P0 = [ Q; anim.t(1:2,:); K0 ]; cnp=6+nK;
  end
  
  if ~isempty(anim.l) % NRSFM
    P0 = [ P0; anim.l ];
    P0 = [ P0(:)' reshape( permute( anim.SBasis, [ 1 3 2 ] ), [], 1 )' ];
    if anim.isProj
      error('has to be an affine camera for NRSFM');
    end
    cnp=6+size(anim.l,1);
    pnp=3*anim.nBasis;
    proj='affineNRSFM';
  else
    P0 = [ P0(:)' anim.S(:)' ];
    if anim.isProj % Projective camera
      if doK
        proj='projectivekap1kap2pp1pp2Ignored';
      else
        proj='projectivek1k2k3k4k5kap1kap2pp1pp2Ignored';
      end
    else % Affine camera
      if doK
        proj='affinetr3k4k5kap1kap2pp1pp2Ignored';
      else
        proj='affinetr3k1k2k3k4k5kap1kap2pp1pp2Ignored';
      end
    end
  end
  fullP = false;
else
  fullP = true;
  if size(anim.P,2)==3
    %homography case
    proj = 'homography';
    P0 = reshape( anim.P, [ 9 anim.nFrame ] );
    P0 = bsxfun(@rdivide, P0(1:8,:), P0(9,:)); cnp=8;
    P0 = [ P0(:)' anim.S(:)' ];
    pnp = 2;
    nFrameFixed = 0;
  else
    if anim.isProj
      P0 = reshape( permute( anim.P, [ 2 1 3 ] ), [ 12 anim.nFrame ] ); cnp=12;
      proj = 'projectiveFull';
    else
      P0 = reshape( permute( anim.P(1:2,:,:), [ 2 1 3 ] ), [ 8 anim.nFrame ] ); cnp=8;
      proj = 'affineFull';
    end
    P0 = [ P0(:)' anim.S(:)' ];
    pnp = 3;
  end
end

% Initialize the point variables
WPermute=permute(anim.W, [1,3,2]);
x=WPermute(repmat(reshape(permute(WMask,[2,1]),1,...
  anim.nFrame,anim.nPoint), [2,1,1])>0);

% launch sba
projac=[ proj 'Jac' projPath ];
proj=[ proj projPath ];

if ~isempty(covx)
  [ ret P info ] = sba( nPoint, 0, nFrame, nFrameFixed, WMask+0.0, ...
    P0, cnp, pnp, x, covx, 2, proj, projac, nItrSBA, 0, [], 'motstr', ...
    K0, KMask );
else
  if ~isempty(anim.l) % NRSFM
    isFirstCoeff1=double(isFirstCoeff1);
    nFrameFixed = 0;
    [ ret P info ] = sba( anim.nPoint, 0, anim.nFrame, nFrameFixed, WMask+0.0, ...
      P0, cnp, pnp, x, 2, proj, projac, nItrSBA, 0, [], 'motstr',isFirstCoeff1,...
      anim.nBasis);
  else
    [ ret P info ] = sba( anim.nPoint, 0, anim.nFrame, nFrameFixed, WMask+0.0, ...
      P0, cnp, pnp, x, 2, proj, projac, nItrSBA, 0, [], 'motstr',...
      KAll, KMask );
  end
end

%  info(1:2)
%info(1:2)/anim.nFrame

% create the output information
if fullP
  if size(anim.P,2)==3
    % homography
    anim.S = reshape( P(8*anim.nFrame + 1 : end ), 2, anim.nPoint );
    P = reshape( P(1:8*anim.nFrame ), [ 8 anim.nFrame ] );
    P(end+1,:) = 1;
    anim.P = reshape( P, [ 3 3 anim.nFrame ] );
  else
    if anim.isProj
      P1 = permute( reshape( P(1:12*anim.nFrame ), [ 4 3 anim.nFrame ] ), [ 2 1 3 ] );
      anim.S = reshape( P(12*anim.nFrame + 1 : end ), 3, anim.nPoint );
    else
      P1 = permute( reshape( P(1:8*anim.nFrame ), [ 4 2 anim.nFrame ] ), [ 2 1 3 ] );
      P1(3,4,:) = 1;
      anim.S = reshape( P(8*anim.nFrame + 1 : end ), 3, anim.nPoint );
    end
    anim.P = P1;
  end
else
  if anim.isProj % Projective camera
    P1 = reshape( P(1:(7+nK)*anim.nFrame ), [], anim.nFrame );
    % the following line is just to prevent automatic update until
    % anim.t, anim.K and anim.R are all full
    anim.t = [];
    anim.K( ~KMask, : ) = P1( 8:end, : );
    anim.R = quaternion( P1(1:4,:) );
    anim.t = P1(5:7,:);
    anim.S = reshape( P(( 7 + nK)*anim.nFrame + 1 : end ), 3, ...
      anim.nPoint );
  else % Affine camera
    if isempty(anim.l)
      P1 = reshape( P(1:(6+nK)*anim.nFrame ), [], anim.nFrame );
      % the following line is just to prevent automatic update until
      % both anim.t and anim.R are full
      anim.R = [];
      anim.t = P1(5:6,:); anim.t(3,:) = 0;
      anim.R = quaternion( P1(1:4,:) );
      if size(P1,1)>=7; anim.K( ~KMask, : ) = P1( 7:end, : ); end
      anim.S = reshape( P(( 6 + nK)*anim.nFrame + 1 : end ), 3, ...
        anim.nPoint );
    else
      P1 = reshape( P(1:cnp*anim.nFrame ), cnp, anim.nFrame );
      % the following line is just to prevent automatic update until
      % both anim.t and anim.R (or anim.l and anim.SBasis) are full
      anim.t = []; anim.SBasis = [];
      anim.R = quaternion( P1(1:4,:) );
      anim.t = P1(5:6,:); anim.t(3,:) = 0;
      anim.l = P1( 7 : end, : );
      anim.SBasis = permute( reshape( P(cnp*anim.nFrame + 1 : end ), 3, ...
        anim.nBasis, anim.nPoint ), [ 1 3 2 ] );
    end
  end
end
