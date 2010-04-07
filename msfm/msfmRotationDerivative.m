function dUdabcd = msfmRotationDerivative(Q)
% wrapper around a mex to find the derivative of a rotation matrix
%
% Find the derivatives of a rotation matrix with respect to its quaternions
% For each quaternion in Q, it returns a [ 2 x 3 x 4 ] matrix: [ 2 x 3 ]
% is the derivative of the first two rows of the rotation, with respect to
% on of the quaternions
%
% USAGE
%  dUdabcd = msfmRotationDerivativeMex(Q);
%
% INPUTS
%  Q          - [ 4 x nFrame ] initial guess of the quaternions
%
% OUTPUTS
%  dUdabcd    - [ 2 x 3 x 4 x nFrame ] optimized quaternions
%
% EXAMPLE
%
% See also
%
% Vincent's Structure From Motion Toolbox      Version 3.0
% Copyright (C) 2009 Vincent Rabaud.  [vrabaud-at-cs.ucsd.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the GPL [see external/gpl.txt]

dUdabcd = msfmRotationDerivativeMex(Q);
