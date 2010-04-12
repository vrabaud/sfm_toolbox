function [ QOpt, err ] = refineExteriorOrientation(S,W,Q)
% wrapper around a mex to optimize the quaternions in exterior orientation
%
% Exterior is : 1 array of 3D position x is known and its 2D orthographic
% projection xp.
% Find the best transformation such that xp=projection*(s*R*x+t)
% Projection is ortho for now
%
% The function here only optimizes the rotation
%
% USAGE
%  Q=refineExteriorOrientation(S,W,Q)
%
% INPUTS
%  S          - [ 3 x nPoint x nFrame ] 3D data
%  W          - [ 2 x nPoint x nFrame ] 2D projection
%  Q          - [ 4 x nFrame ] initial guess of the quaternions
%
% OUTPUTS
%  Q          - [ 4 x nFrame ] optimized quaternions
%  err        - [ 1 x nFrame ] final errors
%
% EXAMPLE
%
% See also
%
% Vincent's Structure From Motion Toolbox      Version 3.0
% Copyright (C) 2009 Vincent Rabaud.  [vrabaud-at-cs.ucsd.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the GPL [see external/gpl.txt]

[ QOpt, err ] = refineExteriorOrientationMex(S,W,Q);
