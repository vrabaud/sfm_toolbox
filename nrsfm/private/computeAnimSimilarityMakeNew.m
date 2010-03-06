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
% Vincent's Structure From Motion Toolbox      Version NEW
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
RRtRR = RR'*RR;
RRtdRRa = RR'*diff(RR,a);
RRtdRRb = RR'*diff(RR,b);
RRtdRRc = RR'*diff(RR,c);
RRtdRRd = RR'*diff(RR,d);

RRtRR = simple(RRtRR);
RRtdRRa = simple(RRtdRRa);
RRtdRRb = simple(RRtdRRb);
RRtdRRc = simple(RRtdRRc);
RRtdRRd = simple(RRtdRRd);

%%%%%  Rest of the code    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
l = [ RRtRR(:) ]
l = [ RRtRR(:), RRtdRRa(:), RRtdRRb(:), RRtdRRc(:), RRtdRRd(:) ];
l = simple( l );

maple restart;
com = 'res := [ ';
for j = 1 : 16
  com = [ com 'Rt[' num2str(j) ']=' char(l(j,1)) ', '];
end
k = 1;
for i = 2 : 5
  for j = 1 : 16
    com = [ com 'dRt[' num2str(k) ']=' char(l(j,i)) ', '];
	k = k + 1;
  end
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

while 1
  tok = regexp( out,'\[(\d*)\]\[(\d*)\]', 'tokens' );
  tokStart = regexp( out,'\[(\d*)\]\[(\d*)\]', 'start' );
  tokEnd = regexp( out,'\[(\d*)\]\[(\d*)\]', 'end' );

  if isempty(tok); break; end

  out = [ out( 1:tokStart(1)-1 ) '(' num2str(str2double(tok{1}{1})+1) ...
      ',' num2str(str2double(tok{1}{2})+1) ')' out( tokEnd(1)+1 : end ) ];
end

fprintf( out );
