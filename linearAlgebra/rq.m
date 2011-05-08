% Perform RQ decomposition
%
% Returns an upper-triangular matrixR and an unitary matrix Q such that
% A=R*Q
%
% REFERENCE HZ2, p579
%
% USAGE
%  [R,Q]=rq(A)
%
% INPUTS
%  A       - input matrix
%
% OUTPUTS
%  R       - upper-triangular matrix
%  Q       - unitary matrix
%
% EXAMPLE
%
% See also
%
% Vincent's Structure From Motion Toolbox      Version 1.1
% Copyright (C) 2008-2011 Vincent Rabaud.  [vrabaud-at-cs.ucsd.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the GPL [see external/gpl.txt]

function [ R Q ] = rq(A)

[Q,R] = qr(A(end:-1:1,end:-1:1)');
Q = Q(end:-1:1,end:-1:1)';
R = R(end:-1:1,end:-1:1)';

if det(Q)<0; R = -R; Q = -Q; end
