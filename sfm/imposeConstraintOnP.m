function [ Q P ] = imposeConstraintOnP( P, isProj, constraint )
% Impose a constraints on Calibrated Projection Matrices
%
% The implemented method is inspired by:
% A Sequential Factorization Method for Recovering Shape and Motion From
% Image Streams (1997)
% Toshihiko Morita, Takeo Kanade
%
% USAGE
%  Q=imposeConstraintOnP( P, isProj )
%
% INPUTS
%  P            - 3x4xT projection matrices
%  isProj       - flag indicating if the camera matrices are projective
%  constraint   - 'ortho' (best Q so that all the rotation parts of P are
%                 as close to ortho as possible), 'orthoHard' (forces the
%                 output rotations of Ps to be ortho), 'firstId' (forces
%                 the first matrix of P to be eye(3,4) ),'firstIdRot'
%                 (forces P(:,1:3,1) to be eye(3) )
%
% OUTPUTS
%  Q            - matrix such that P(:,:,i)*Q contains a rotation matrix
%                 Also, the first matrix is the identity
%  P            - new projection matrix, with Q applied, such that the
%                 first matrix is identity
%
%
% EXAMPLE
%
% See also
%
% Vincent's Structure From Motion Toolbox      Version 1.1
% Copyright (C) 2009 Vincent Rabaud.  [vrabaud-at-cs.ucsd.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the GPL [see external/gpl.txt]

nFrame=size( P, 3 );

% Compute the orthonormality constraints
switch constraint
  case 'ortho'
    if isProj
      G=zeros( 6*nFrame, 6 );
      for i=1:nFrame
        G(i,:) = g( P(1,:,i), P(1,:,i) );
        G(nFrame+i,:) = g( P(2,:,i), P(2,:,i) );
        G(2*nFrame+i,:) = g( P(2,:,i), P(2,:,i) );
        G(3*nFrame+i,:) = g( P(1,:,i), P(2,:,i) );
        G(4*nFrame+i,:) = g( P(1,:,i), P(3,:,i) );
        G(5*nFrame+i,:) = g( P(2,:,i), P(3,:,i) );
      end
      c = [ ones(3*nFrame,1); zeros(3*nFrame,1) ];
    else
      G=zeros( 3*nFrame, 6 );
      for i=1:nFrame
        G(i,:) = g( P(1,:,i), P(1,:,i) );
        G(nFrame+i,:) = g( P(2,:,i), P(2,:,i) );
        G(2*nFrame+i,:) = g( P(1,:,i), P(2,:,i) );
      end
      c = [ ones(2*nFrame,1); zeros(nFrame,1) ];
    end
    
    % Solve for the square root matrix
    l = G\c;
    
    L = [ l(1:3)'; l(2) l(4:5)'; l(3) l(5:6)' ];
    
    [ U S V ] = svd(L);
    S( S<0 ) = 0;
    
    Q = U*sqrt(S)*V';
    Q(4,4) = 1;
  case 'orthoHard'
    [ Q P ] = imposeConstraintOnP( P, isProj, 'ortho' );
  case 'firstId'
    % Compute the new projection matrices
    Q = imposeConstraintOnP( P, isProj, 'firstIdRot' );
    
    if isProj
      Q(1:3,4) = P(:,1:3,1)\(-P(:,4,1));
    else
      Q(1:3,4) = pinv( P(:,1:3,1) )*([ 0; 0; 1] - P(:,4,1));
    end
  case 'firstIdRot'
    % Compute the new projection matrices
    Q = zeros(3,3);
    if isProj
      temp = [ 2 3; 1 3; 1 2 ];
      for i=1:3
        Q(:,i) = vect( cross( P(temp(i,1),1:3,1)', ...
          P(temp(i,2),1:3,1) ), 'h');
        Q(:,i) = Q(:,i)/(P(i,1:3)*Q(1:3,i));
      end
    else
      temp = [ 2 1 ];
      Q(:,3) = cross( P(1,1:3,1), P(2,1:3,1) )';
      for i=1:2
        Q(:,i) = vect( cross( P(temp(i),1:3,1)', Q(:,3) ), 'h' );
        Q(:,i) = Q(:,i)/(P(i,1:3)*Q(:,i));
      end
    end
    if det(Q)<0; Q(:,3) = -Q(:,3); end
    Q(4,4) = 1;
end

if nargout==2
  switch constraint
    case 'ortho',
      P(:,1:3,:)=multiTimes(P(:,1:3,:),Q(1:3,1:3),1);
    case 'orthoHard',
      if isProj
        for i=1:nFrame
          P(:,1:3,i) = rotationMatrix( P(:,1:3,i) );
        end
      else
        for i=1:nFrame
          temp = rotationMatrix( rotationMatrix( P(1:2,1:3,i ) ) );
          P(1:2,1:3,i) = temp(1:2,:);
        end
      end
    case 'firstId',
      P=multiTimes(P,Q,1);
      P(:,:,1) = eye(3,4); if ~isProj; P(3,3:4,1) = [0 1]; end
    case 'firstIdRot',
      P(:,1:3,:)=multiTimes(P(:,1:3,:),Q(1:3,1:3),1);
      P(:,1:3,1) = eye(3,3); if ~isProj; P(3,3,1) = 0; end
  end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function res=g(a,b)
res=[ a(1)*b(1) a(1)*b(2)+a(2)*b(1) a(1)*b(3)+a(3)*b(1) a(2)*b(2) ...
  a(2)*b(3)+a(3)*b(2) a(3)*b(3) ];
end
