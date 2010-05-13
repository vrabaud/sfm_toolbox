% isInAData contains 1 or 0 if the following parameters are in adata:
%  q1, q2, q3, q4, tr1, tr2, tr3, k1, k2, k3, k4, k5, kap1, kap2, pp1, pp2
% camType: 0 is the full camera matrix, 1 is K,R,t, 2 is a homography
function str = makeSBASlave( fid, camType, isProj, isCalibrated, hasDistortion )

reset(symengine);

switch camType
  case 0,
    varCam = { 'p1' 'p2' 'p3' 'p4' 'p5' 'p6' 'p7' 'p8' 'p9' 'p10' 'p11' ...
      'p12' };
    varPoint = { 'x' 'y' 'z' };
    syms x y z real;
    syms p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 real;

    if isProj
      fileName = 'projectiveFull';
    else
      fileName = 'affineFull';
      p9 = 0; p10 = 0; p11 = 0; p12 = 1;
    end

    P = [ p1 p2 p3 p4; p5 p6 p7 p8; p9 p10 p11 p12 ];
    proj = P*[ x y z 1 ]';
    proj = proj(1:2)/proj(3);
    camParamStr = 'p'; pointParamStr = 'xyz';
  case 2,
    varCam = { 'h1' 'h2' 'h3' 'h4' 'h5' 'h6' 'h7' 'h8' };
    varPoint = { 'x' 'y' };
    syms x y real;
    syms h1 h2 h3 h4 h5 h6 h7 h8 real;

    H = [ h1 h4 h7; h2 h5 h8; h3 h6 1 ];
    proj = H*[ x y 1 ]';
    proj = proj(1:2)/proj(3);
    fileName = 'homography';
    camParamStr = 'h'; pointParamStr = 'xyz';
  case 1,
    varCam = { 'q1' 'q2' 'q3' 'q4' 'tr1' 'tr2' 'tr3' 'k1' 'k2' 'k3' 'k4'...
      'k5' 'kap1' 'kap2' 'pp1' 'pp2' };
    varPoint = { 'x' 'y' 'z' };
    syms x y z real;
    syms q1 q2 q3 q4 tr1 tr2 tr3 k1 k2 k3 k4 k5 kap1 kap2 pp1 pp2 real;

    % Some symbolic math
    Pi0 = eye(3,4); 
    if isProj
      fileName = 'projective';
    else
      Pi0(3,3:4) = [ 0 1 ]; tr3 = 0; k4 = 0; k5 = 0;
      fileName = 'affine';
    end
    
    pointParamStr = 'xyz';
    if isCalibrated
      k1 =1; k2 = 0; k3 = 1; k4 = 0; k5 = 0; camParamStr = 'rt';
    else
      camParamStr = 'rtk';
    end
    if hasDistortion
      camParamStr = [ camParamStr 'dis' ];
    else
      kap1 = 0; kap2 =0; pp1 =0; pp2 = 0;
    end

    % Compute the rotation matrix from the quaternions
    K = [ k1 k2 k4; 0 k3 k5; 0 0 1 ];
    R = [ q1^2+q2^2-q3^2-q4^2 2*q2*q3-2*q1*q4 2*q1*q3+2*q2*q4; ...
      2*q1*q4+2*q2*q3 q1^2-q2^2+q3^2-q4^2 2*q3*q4-2*q1*q2; ...
      2*q2*q4-2*q1*q3 2*q1*q2+2*q3*q4 q1^2-q2^2-q3^2+q4^2 ]/...
      (q1^2+q2^2+q3^2+q4^2);

    proj = Pi0*[ R [tr1 tr2 tr3]'; 0 0 0 1]*[ x y z 1 ]';
    proj = proj(1:2)/proj(3);

    % Apply distortion
    % [x]     [x]
    % [y] = R*[y] + t
    % [z]     [z]
    %
    % x' = x/z
    % y' = y/z
    %
    % x" = x'*(1 + k,,1,,r^2^ + k,,2,,r^4^) + 2*p,,1,,x'*y' + p,,2,,(r^2^+2*x'^2^)
    % y" = y'*(1 + k,,1,,r^2^ + k,,2,,r^4^) + p,,1,,(r^2^+2*y'^2^) + 2*p,,2,,*x'*y'
    % where r^2^ = x'^2^+y'^2^
    %
    % u = fx*x" + cx
    % v = fy*y" + cy

    r2 = sum( proj.^2 );
    r4 = r2^2;
    proj = [ proj(1)*(1 + kap1*r2 + kap2*r4) + 2*pp1*proj(1)*proj(2) + ...
      pp2*(r2+2*proj(1)^2), ...
      proj(2)*(1 + kap1*r2 + kap2*r4) + pp1*(r2+2*proj(2)^2) + ...
      2*pp2*proj(1)*proj(2), 1 ];
    proj = K(1:2,:)*proj';
end

% Simplify the whole projection expression
proj = simpleArray( subs(proj') );

% Compute the Jacobian expressions
sym A real;
for i = 1 : length(varPoint)
  try eval( [ 'tmp = diff(proj,' varPoint{i} ');' ] );
  catch; tmp = 0;
  end
  A(:,i) = tmp;
end

isPresent = zeros(1, length(varCam));
for i = 1 : length(varCam)
  if ~eval( [ 'isfloat(' varCam{i} ')' ] )
    isPresent(i) = 1;
    try eval( [ 'tmp = diff(proj,''' varCam{i} ''');' ] );
    catch; tmp = 0;
    end
    A(:,end+1) = tmp;
  end
end

% A = simpleArray(A); % Actually reduces the codegen optimization

% Get the name of the function
if camType==1
  if any(isPresent==0)
    fileName = [ fileName varCam{ isPresent==0 } 'Ignored' ];
  end
end

%%%%% Projection
% j is camera, i point !
str{1} = [ 'API_MOD void CALL_CONV ' fileName '(int j, double *' ...
  camParamStr ', double *' pointParamStr ', double *xij, double **adata) { \n' ];

str{2} = 'double ';
for i = 1 : length(varPoint); str{2} = [ str{2} varPoint{i} ', ' ]; end
for i = 1 : length(varCam)
  if isPresent(i); str{2} = [ str{2} varCam{i} ', ' ]; end
end
str{2} = [ str{2}(1:end-2) ';\n\n' ];

% Get adata parameters if any
str{3} = '';
if camType==1
  ind = 0;
  if ~isCalibrated || hasDistortion
    str{3} = [ str{3} 'int ind = ' int2str(6+isProj) ';\n\n' ];
  end
  if ~isCalibrated
    str{3} = [ str{3} 'double *k = adata[0];\n' ];
    str{3} = [ str{3} 'double *kMask = adata[1];\n' ];
    ind = 2;
  end
  if hasDistortion
    str{3} = [ str{3} 'double *distor = adata[' int2str(ind) '];\n' ];
    str{3} = [ str{3} 'double *distorMask = adata[' int2str(ind+1) '];\n' ];
  end
end

% Define variables
j = 0;
% Get the point parameters
for i = 1 : length(varPoint)
    str{3} = [ str{3} varPoint{i} ' = ' pointParamStr '[' int2str(i-1) '];\n' ];
end
% Get the camera parameters
for i = 1 : length(varCam)
  if camType==0 || camType==2 || (camType==1 && i<=7)
    if isPresent(i)
      str{3} = [ str{3} varCam{i} ' = ' camParamStr '[' int2str(i-1) '];\n' ];
    end
  else % read from adata
    ii = i-8;
    if isPresent(i) 
      if ii<=4 % K param
        str{3} = [ str{3} varCam{i} ' = kMask[' int2str(ii) '] ? k[5*j+' int2str(ii) '] : ' camParamStr '[ind++];\n' ];
      else %distor param
        str{3} = [ str{3} varCam{i} ' = distorMask[' int2str(ii-5) '] ? distor[4*j+' int2str(ii-5) '] : ' camParamStr '[ind++];\n' ];
      end
    end
  end
end
str{3} = [ str{3} '\n' ];

% Compute the projection
str{4} = [ toC( proj, 'xij' ) '\n}\n\n' ];


%%%%% Jacobian
str{4} = [ str{4} 'API_MOD void CALL_CONV ' fileName 'Jac(int j, double *' ...
  camParamStr ', double *' pointParamStr ', double *Aij, double *Bij, double **adata) { \n' ];

str{5} = str{2};
str{6} = str{3};


j = 2*length(varPoint);
A = [ reshape(A(:,1:j/2)',[],1)' reshape(A(:,j/2+1:end)',[],1)' ];
str{7} = [ toC( reshape(A,1,[]), 'Aij' ) ];

for i = 0 : j-1
  str{7} = strrep( str{7}, [ 'Aij[' int2str(i) ']' ], [ 'Bij[' int2str(i) ']' ] );
end
for i = j : length(A)
  str{7} = strrep( str{7}, [ 'Aij[' int2str(i) ']' ], [ 'Aij[' int2str(i-j) ']' ] );
end

% Only fill the existing entries in Aij
varExtraCamNbr = sum(isPresent(8:end));
if camType==1 && varExtraCamNbr>0
  str{7} = [ 'ind = ' int2str(6+isProj) ';\n' str{7} ];
  for i = 8 : length(varCam)
    ii = i-8;
    if isPresent(i) 
      if ii<=4 % K param
        mask = 'kMask';
      else %distor param
        ii = ii - 5; mask = 'distorMask';
      end
      for k = i-1 + [ 0 6+isProj+varExtraCamNbr ] - (~isProj)
        str{7} = regexprep( str{7}, [ 'Aij\[' int2str(k) '\] = (.*?);' ], [ 'if (!' mask '[' int2str(ii) ']) \n    Aij[ind++] = $1;' ] );
      end
    end
  end
  for i = 6+isProj+varExtraCamNbr + [ 0 : 7 ]
    str{7} = regexprep( str{7}, [ 'Aij\[' int2str(i) '\] = (.*?);' ], 'Aij[ind++] = $1;\n' );
  end
end

str{7} = [ str{7} '}\n\n' ];

% Deal with the temporary variables (t80 ...)
for n=[4 7]
  tmp = regexp( str{n}, 't(\d+)', 'tokens' );

  strt = 'double ';
  isDone = zeros(1,1000);
  for i=1:length(tmp)
    tmp2 = round( str2double(tmp{i}{1}) );
    if tmp2>length(isDone) || ~isDone(tmp2)
      strt = [ strt 't' tmp{i}{1} ', ' ];
      isDone( tmp2 ) = 1;
    end
  end
  if ~isempty(tmp)
    str{n-2} = [ str{n-2} strt(1:end-2) ';\n\n' ];
  end
end

% Write everything to a file
for i=1:length(str)
  fprintf(fid,str{i});
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
function out = toC( a, var )
  str = char(collect(a));
  out = evalin(symengine, [ 'opts := generate::optimize(expr(' str '));' ] );
  evalin(symengine, [ var ':=rhs(opts[-1]);' ] );

  out = [ char(evalin(symengine, [ 'generate::C(opts[1..-2]);' ] )) ...
    char( evalin(symengine, [ 'generate::C(' var ')' ] )) ];
  out = strrep( out, '"', '' );
  out = strrep( out, [ var '[0][' ], [ var '[' ] );
end

function a = simpleArray( a )
  for jSub = 1 : size(a,1)
    for iSub = 1 : size( a, 2 )
      try
        a(jSub,iSub) = simple( collect( a(jSub,iSub) ) );
      catch
      end
    end
  end
end
