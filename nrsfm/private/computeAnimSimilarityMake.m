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

%%%%%  Infimum in TSFM   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms a b c cc d dd real

R = [ a^2+b^2-c^2-d^2 2*b*c-2*a*d 2*a*c+2*b*d; ...
 2*a*d+2*b*c a^2-b^2+c^2-d^2 2*c*d-2*a*b; ...
 2*b*d-2*a*c 2*a*b+2*c*d a^2-b^2-c^2+d^2 ]/(a^2+b^2+c^2+d^2);
R = R(1:2,:);
R1 = simple( R(:,3)*R(:,3)'/(R(:,3)'*R(:,3)) - eye(2) );
R2 = simple( -R1*R(:,1:2) );
RR = subs( [ R1 R2 ] );

%%%%%  Code for optimal 2-view SFM  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms a b c cc d dd real

R = [ a^2+b^2-c^2-d^2 2*b*c-2*a*d 2*a*c+2*b*d; ...
 2*a*d+2*b*c a^2-b^2+c^2-d^2 2*c*d-2*a*b; ...
 2*b*d-2*a*c 2*a*b+2*c*d a^2-b^2-c^2+d^2 ]/(a^2+b^2+c^2+d^2);
Rbar = R;
% Rbar(:,3) = -R(3,:);              Used to show the equivalence
% Rbar(1:2,1:2)=R(1:2,1:2)';

Rbar = Rbar(1:2,:);
Pi = eye(2,3);
RbarInv = ( Pi'*Pi + Rbar'*Rbar )^(-1);
Rbar1 = simple( eye(2) - Pi*RbarInv*Pi' );
Rbar2 = simple( -Pi*RbarInv*Rbar' );
RR2 = subs( [ Rbar1 Rbar2 ] );

% 1/2*RR1+RR2

%%%%%  Code for optimal 3-view SFM  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  syms a b c d aa bb cc dd real
%  clear R;
%  R{1} = eye(2,3);
%  R{2} = [ a^2+b^2-c^2-d^2 2*b*c-2*a*d 2*a*c+2*b*d; ...
%   2*a*d+2*b*c a^2-b^2+c^2-d^2 2*c*d-2*a*b; ...
%   2*b*d-2*a*c 2*a*b+2*c*d a^2-b^2-c^2+d^2 ]/(a^2+b^2+c^2+d^2);
%  R{3} = [ aa^2+bb^2-cc^2-dd^2 2*bb*cc-2*aa*dd 2*aa*cc+2*bb*dd; ...
%   2*aa*dd+2*bb*cc aa^2-bb^2+cc^2-dd^2 2*cc*dd-2*aa*bb; ...
%   2*bb*dd-2*aa*cc 2*aa*bb+2*cc*dd aa^2-bb^2-cc^2+dd^2 ]/(aa^2+bb^2+cc^2+dd^2);
%  for i=2:3; R{i}=R{i}(1:2,:); end
%  RInv = 0; for i=1:3; RInv = RInv + R{i}'*R{i}; end; RInv=simple(inv(RInv)); 
%  RR = sym(zeros( 2, 6, 3 ));
%  
%  for i=1:3; for j=1:3; RR(:,2*j-1:2*j,i)=(i==j)*eye(2)-R{i}*RInv*R{j}'; end; end
% Way too complex expressions

%%%%%  Rest of the code    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
l = RR;
l = [ l; diff(l,a); diff(l,b); diff(l,c); diff(l,d) ];
l = [ l; diff(l(3:end,:),a); diff(l(5:end,:),b); diff(l(7:end,:),c); ...
    diff(l(9:end,:),d) ];
l = simple( l' );

% maple restart;
% com = [ 'res := [ val=' char(l(1)) ', ' ];
% for j = 1 : 4
%   com = [ com 'dRes[' num2str(j) ']=' char(l(j+1)) ', '];
% end
% com = [ com(1:end-2) '];' ];

%  %%% Original version of optimization: ok but there is better
%  in = [ 'l:=' char( l ) ];
%  
%  maple( in );
%  out = maple( 'lOpt := codegen[C](l,optimized)' );

%%% Now, use the following:
%%% You need an expression of the form : '[ a[1][1] = x, a[1][2] =y^2 ]'
%%% not the overall brackets and the indexing starting from 1 (as it is
%%% maple. The quotes: simply coz it has to be a string
%%% so you create an array of assigments basically
maple restart;
com = 'res := [ ';
for i = 1 : 2
  for j = 1 : 4
    com = [ com 'R[' num2str(j) '][' num2str(i) ']=' char(l(j,i)) ', '];
  end
end
for i = 3 : 10
  for j = 1 : 4
    com = [ com 'dR[' num2str(j) '][' num2str(i-2) ']=' char(l(j,i)) ', '];
  end
end
for i = 11:30
  for j = 1 : 4
    com = [ com 'ddR[' num2str(j) '][' num2str(i-10) ']=' char(l(j,i)) ', '];
  end
end
com = [ com(1:end-2) '];' ];

%%%%%  Code for optimal Rotation in Pose%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms a b c d w1 w2 x y z real

R = [ a^2+b^2-c^2-d^2 2*b*c-2*a*d 2*a*c+2*b*d; ...
 2*a*d+2*b*c a^2-b^2+c^2-d^2 2*c*d-2*a*b; ...
 2*b*d-2*a*c 2*a*b+2*c*d a^2-b^2-c^2+d^2 ]/(a^2+b^2+c^2+d^2);
R = R(1:2,:);

l=[ w1; w2 ]- R*[x;y;z]; l=l'*l;
l = [ l; diff(l,a); diff(l,b); diff(l,c); diff(l,d) ];
l = simple( l' );

maple restart;
com = [ 'res := [ val=' char(l(1)) ', ' ];
for j = 1 : 4
  com = [ com 'dVal[' num2str(j) ']=' char(l(j+1)) ', '];
end
com = [ com(1:end-2) '];' ];

%%%%%  Below is the code to generate C-code    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

while 1
  tok = regexp( out,'\[(\d*)\]\[(\d*)\]', 'tokens' );
  tokStart = regexp( out,'\[(\d*)\]\[(\d*)\]', 'start' );
  tokEnd = regexp( out,'\[(\d*)\]\[(\d*)\]', 'end' );

  if isempty(tok); break; end

  out = [ out( 1:tokStart(1)-1 ) '(' num2str(str2double(tok{1}{1})+1) ...
      ',' num2str(str2double(tok{1}{2})+1) ')' out( tokEnd(1)+1 : end ) ];
end

fprintf( out );
