function A = skew(a)
% Return the skew matrix or its original vector
%
% USAGE
%  A=skew(a)
%
% INPUTS 1 - returns the skew matrix of a vector
%  a       - 3x1 or 1x3 vector
%
% INPUTS 2 - returns the vector that created the closest skew matrix
%  a       - 3x3 skew matrix
%
% OUTPUTS 1
%  A       - corresponding skew matrix
%
% OUTPUTS 2
%  A       - vector that created the matrix
%
% EXAMPLE
%
% See also
%
% Vincent's Structure From Motion Toolbox      Version 3.1.1
% Copyright (C) 2008-2011 Vincent Rabaud.  [vrabaud-at-cs.ucsd.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the GPL [see external/gpl.txt]

if numel(a)==3
  % returns the skew matrix of a vector
  % Reference: HZ2, p581, equation (A4.5)
  A=[0 -a(3) a(2); a(3) 0 -a(1); -a(2) a(1) 0];
  return
end

if all(size(a)==[3,3])
  % returns the vector that created the closest skew matrix
  A=0.5*[ a(3,2)-a(2,3); a(1,3)-a(3,1); a(2,1)-a(1,2)];
  return
end

error('Bad input dimensions. Input must be 1x3, 3x1 or 3x3');
