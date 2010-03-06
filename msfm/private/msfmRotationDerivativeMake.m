% Function that optimizes Pose given an initial solution of the rotation only
%
% IMPORTANT: THIS FUNCTION ONLY DEALS WITH THE ORTHOGRAPHIC CASE RIGHT NOW
%
% This function optimizes the pose estimation of one camera with respect to
% another given only the 2 corresponding views.  The data is supposed to be
% centered and only the rotation has to be estimated. I it is here 
% performed using simple gradient descent.
%
% The following code explains how the code to compute the
% error/gradient/hessian was created using matlab symbolic toolbox.
%
% Vincent's Structure From Motion Toolbox      Version 2.0\n
% Copyright (C) 2009 Vincent Rabaud.  [vrabaud-at-cs.ucsd.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the Lesser GPL [see external/lgpl.txt]

%%%%%  Define the rotation matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms a b c d real

R = [ a^2+b^2-c^2-d^2 2*b*c-2*a*d 2*a*c+2*b*d; ...
 2*a*d+2*b*c a^2-b^2+c^2-d^2 2*c*d-2*a*b; ...
 2*b*d-2*a*c 2*a*b+2*c*d a^2-b^2-c^2+d^2 ]/(a^2+b^2+c^2+d^2);
R = R(1:2,:);

%%%%%  Compute its derivatives with respect to the quaternions %%%%%%%%%%%%%%%%%%%%%%
l = R;
l = [ diff(l,a) diff(l,b) diff(l,c) diff(l,d) ];
l = simple( l );

%%% Now, use the following:
%%% You need an expression of the form : '[ a[1][1] = x, a[1][2] =y^2 ]'
%%% not the overall brackets and the indexing starting from 1 (as it is
%%% maple. The quotes: simply coz it has to be a string
%%% so you create an array of assigments basically
maple restart;
com = 'res := [ ';
k = 1;
for i = 1 : 12
  for j = 1 : 2
    com = [ com 'dR[' num2str(k) ']=' char(l(j,i)) ', '];
	k = k + 1;
  end
end
com = [ com(1:end-2) '];' ];

%%%%%  Below is the code to generate C-code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Mupad try, but it sucks ...
%  reset(symengine);
%  out = evalin(symengine, [ 'opts := generate::optimize([R=expr(' char(Rp') ') ]);' ] );
%  
%  out = evalin(symengine, [ 'opts := generate::optimize([R=expr(' char(Rp') '), dR=expr(' char(l(:,3:10)) '), ddR=expr(' char(l(:,11:end)) ')]);' ] );
%  evalin(symengine, [ ':=rhs(opts[-1]);' ] );
%  
%  
%  
%  out = [ char(evalin(symengine, [ 'generate::C(opts[1..-2]);' ] )) ...
%    char( evalin(symengine, [ 'generate::C(rhs(opts[-1]))' ] )) ];
%  out = strrep( out, '"', '' );
%  out = strrep( out, [ var '[0][' ], [ var '[' ] );


%%% now, save your crazy array into the res variable
maple( com );

%  %%% cost of the assigments with no optimization
%  maple('codegen[cost](res)')
%  % maple( 'opt1 := [codegen[optimize](l)]' );
%  %%% cost of the assigment with the normal optimization
%  in = [ 'l:=' char( l ) ]; maple( in ); maple('codegen[cost](codegen[optimize](l))')

%%% code for better optimization and corresponding cost
maple( 'opt2 := [codegen[optimize](res,tryhard)]' );
maple('codegen[cost](opt2)')

%%% Convert to C code
out = maple( 'codegen[C](opt2)' );

%%% clean the C code and convert it to Matlab code (replace [] by ()
out = strrep( out, ';   ', ';\n' );
out = strrep( out, '~', '' );


tout = ''
existing = cell(1,0);
tmp = regexp( out,'(t\d*)', 'tokens' )
for i = 1 : length(tmp)
  ii = tmp(i);
  ii = ii{1}{1};
  doExist = false;
  for j = 1 : size(existing,2)
    if strcmp(existing{j},ii)
	  doExist = true;
	  break;
	end
  end
  if ~doExist
    existing{1,end+1} = ii;
    tout = [ tout, ', ', ii ];
  end
end
out = [ 'double ', tout(2:end), ';\n', out ];
fprintf(out)
