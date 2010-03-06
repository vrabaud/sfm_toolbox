% Transform a matrix into a horizontal or vertical vector
%
% USAGE
%  x = vect( x, ori )
%
% INPUTS
%  x       - vector or matrix
%  ori     - orientation ( 'v' or 'h' )
%
% OUTPUTS
%  x       - vectorized x in the right orientation
%
% EXAMPLE
%  x = zeros( 1, 0 );
%  x(:)
%  vect( zeros(0,1), 'v' )
%
% See also
%
% Vincent's Structure From Motion Toolbox      Version 1.1
% Copyright (C) 2009 Vincent Rabaud.  [vrabaud-at-cs.ucsd.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the Lesser GPL [see external/lgpl.txt]

function x = vect( x, ori )

if strcmp(ori,'h'); x=x'; end
x = x(:);

if strcmp(ori,'v') % vertical, n by 1
  if size(x,1)==1; x = x'; end
else % horizontal, 1 by n
  if size(x,2)==1; x = x'; end
end
