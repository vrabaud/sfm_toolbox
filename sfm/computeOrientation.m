function [ R t s ]=computeOrientation(x,xp,method,varargin)
% Compute Absolute or Exterior Orientation (Pose Estimation)
%
% Absolute is : given x and xp (2 arrays of 3D coordinates), find the best
% transformation such that x=s*R*xp+t
%
% Exterior is : 1 array of 3D position x is known and its 2D orthographic
% projection xp.
% Find the best transformation such that xp=projection*(s*R*x+t)
% (same as Pose Estimation, ePNP). Projection is ortho for now
%
% The routines below are only for the orthographic case for now
%
% USAGE
%  [R,t,s]=computeOrientation(x,xp,method)
%
% INPUTS
%  x,xp       - 3xN or 2xN array of points or right/left.
%               For 'exteriorSequence', it can be of size
%               2 x nPoint x nFrame and 3 x nPoint x nFrame
%  method     - 'absolute', 'absoluteHard' (try more if Horn's fails)
%               or 'exterior' and 'exteriorSequence'
%  varargin   - list of paramaters in quotes alternating with their values
%       - 'RIni' initial rotation from which the optimization will be
%       started (only for 'exterior' and 'exteriorSequence')
%
% OUTPUTS
%  R       - rotation matrix
%  t       - translation vector
%  s       - scale factor (if not requested, s=1)
%
% EXAMPLE
%
% See also
%
% Vincent's Structure From Motion Toolbox      Version 3.1
% Copyright (C) 2008-2011 Vincent Rabaud.  [vrabaud-at-cs.ucsd.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the GPL [see external/gpl.txt]

nPoint=size(x,2);

switch method
  case { 'absolute', 'absoluteHard' }
    % Computes the absolute orientation between 2 sets of 3D points
    % Xright=s*R*Xleft + t   . We want to recover s,R,T
    % Ref: B.K.P. Horn, H.M. Hilden, and S. Negahdaripour, Closed-Form
    % Solution of Absolute Orientation Using Orthonormal Matrices
    rr=normalizePoint(x,4); rl=normalizePoint(xp,4);
    rrBar=mean(rr,2); rlBar=mean(rl,2);
    rrp=rr-rrBar(:,ones(1,nPoint)); rlp=rl-rlBar(:,ones(1,nPoint));
    
    M=zeros(3); for i=1:nPoint; M=M+rrp(:,i)*rlp(:,i)'; end
    [V,D]=eig(M'*M);
    V=real(V);
    
    R=1/sqrt(D(2,2))*V(:,2)*V(:,2)' + 1/sqrt(D(3,3))*V(:,3)*V(:,3)';
    temp=V(:,1)*V(:,1)';
    if D(1,1)>0
      R=R+1/sqrt(D(1,1))*temp;
    else
      if det(R+temp)>0; R=R+temp; else R=R-temp; end
    end
    R=M*R;
    
    R=rotationMatrix(R);
    
    if det(R)<0 && strcmp(method,'absoluteHard') % if Horn's method fails
      warning('Horn''s method failed. Trying GloptiPoly3 solution');
      RHorn=R;
      try
        % Sedumi and GloptiPoly3 must be installed and in the path !
        mpol R 3 3;
        
        g0=rlp-R*rrp;
        g0=g0(1,:)*g0(1,:)'+g0(2,:)*g0(2,:)'+g0(3,:)*g0(3,:)';
        
        % define the rotation constraints
        K = [ R(1,:)*R(1,:)' == 1, R(2,:)*R(2,:)'==1, R(3,:)*R(3,:)'==1,...
          R(1,:)*R(2,:)' == 0, R(1,:)*R(3,:)' == 0, R(2,:)*R(3,:)' == 0,...
          det(R)>=0 ];
        
        % define the problem and solve it
        P=msdp(min(g0),K);
        [ status obj ] = msol(P);
        
        R=rotationMatrix( double(R) );
      catch %#ok<CTCH>
        R = RHorn;
      end
    end
    
    % Figure out s and t
    if nargout==2
      s=1;
    else
      s=norm(rrp,'fro')/norm(rlp,'fro');
      if norm( rrp+s*R*rlp, 'fro' ) < norm( rrp-s*R*rlp, 'fro' ); s=-s; end
    end
    
    t=rrBar-s*R*rlBar;
  case 'exterior'
    if ~isempty(varargin); RIni = varargin{0};
    else RIni=[];
    end

    if isempty(RIni)
      try
        % Sedumi and GloptiPoly3 must be installed and in the path !
        mpol R 2 3;
        if nargout<2; t=0; else mpol t 2 1; end
        if nargout<3; s=1; else mpol s; end
        
        g0=0;
        for i=1:nPoint
          tmp=xp(:,i) - s*(R*x(:,i)+t); %#ok<NODEF>
          g0 = g0 + tmp'*tmp;
        end
        
        % define the rotation constraints
        K = [ R(1,:)*R(1,:)' == 1, R(2,:)*R(2,:)'==1, R(1,:)*R(2,:)' == 0];
        
        % define the problem and solve it
        mset('verbose',false);
        P=msdp(min(g0),K);

        [ status obj ] = msol(P);
        
        R=rotationMatrix( double(R) );
        
        if nargout>=2; t=[ double(t); 0]; end
        if nargout==3; s=double(s); end
      catch %#ok<CTCH>
        warning(['GloptiPoly3 not installed or Gloptypoly crashed,' ...
          'using ePnP']);
        try
          mset clear;
        catch %#ok<CTCH>
          [ R t ] = efficient_pnp( x', xp', eye(3,3) );
          t=t(1:2);
        end
      end
      if any(isnan(R))
        warning(['GloptiPoly failed, using ePnP']);
        [ R t ] = efficient_pnp( x', xp', eye(3,3) );
        t=t(1:2);
      end
    else
      R=RIni;
    end
    
    % perform gradient descent to optimize the rotation
    Q = refineExteriorOrientation(x,xp,quaternion(R));
    R=quaternion( Q );
  case 'exteriorSequence'
    RIni=getPrmDflt( varargin, {'RIni' [] }, 1 );
    if length(RIni)>1
      Q=quaternion(RIni);
    else
      Q=2*rand(4,size(xp,3))-1;
    end
    [ Q err ]=refineExteriorOrientation(x,xp,Q);
    
    for i = 1 : 10
      [ QNew errNew ]=refineExteriorOrientation(x,xp,2*rand(4,size(xp,3))-1);
      temp = errNew<err;
      Q(:,temp) = QNew(:,temp);
    end
    
    R=quaternion( Q );
end
end
