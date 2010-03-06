function [ Sim QTot ] = computeCsfmInfimum( W )
% Compute theinfimum for all pairs of frames Bundle Adjustment
%
% Just a wrapper aroundthe mex file
%
% USAGE
%  [ Sim QTot ] = computeCsfmInfimum( W )
%
% INPUTS
%  W      - [ 2 x nPoint x nFrame ] measurements
%
% OUTPUTS
%  Sim    - [ nFrame x nFrame ] infimum matrix for every pair of frames
%  QTot   - [ 4 x nFrame x nFrame ] optimal quaternions
%
% EXAMPLE
%
% See also
%
% Vincent's Structure From Motion Toolbox      Version NEW
% Copyright (C) 2009 Vincent Rabaud.  [vrabaud-at-cs.ucsd.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the Lesser GPL [see external/lgpl.txt]

[ Sim QTot ] = computeCsfmInfimumMex( W );

