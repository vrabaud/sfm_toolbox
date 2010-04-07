function [ H nMax ] = computeHomography( x1, x2, method, thres, Sim )
% Compute homography between two sets of 2D points
%
% The computed homography is such that x2=H*x1
%
% USAGE
%  H = computeHomography(x1,x2,1)
%  H = computeHomography(x1,x2,10,thres)
%  H = computeHomography(x1,x2,Inf,[],Sim)
%
% INPUTS
%  x1       - [ 2 x nPoint ] or [ 3 x nPoint ] first set of 2D points
%  x2       - [ 2 x nPoint ] or [ 3 x nPoint ] second set of 2D points
%  method   - 1 for 4-point algorithm, 10 for RANSAC and Inf for PROSAC
%  thres    - matching threshold for RANSAC
%  Sim      - similarity matrix for PROSAC
%
% OUTPUTS
%  H        - best homography such that x2 = H * x1
%
% EXAMPLE
%
% See also
%
% Vincent's Structure From Motion Toolbox      Version 2.1
% Copyright (C) 2009 Vincent Rabaud.  [vrabaud-at-cs.ucsd.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the GPL [see external/gpl.txt]

%  x1=normalizePoint(x1,3); x2=normalizePoint(x2,3);

switch method
  case 1
    % Reference: HZ2, p.109, Algorithm 4.2
    nPoint=size(x1,2);

    [ xT1 T1 ]=normalizePoint(x1,-Inf);
    [ xT2 T2 ]=normalizePoint(x2,-Inf);

    % (i) (ii)
    switch size(x1,1)
      case 2, % common DLT
        A=zeros(2*nPoint,9);
        A(1:2:end,4:6) = -xT1';
        A(1:2:end,7:9) = (repmat(xT2(2,:),[3 1]).*xT1)';
        A(2:2:end,1:3) = xT1';
        A(2:2:end,7:9) = -(repmat(xT2(1,:),[3 1]).*xT1)';
      case 3,
        A=zeros(3*nPoint,16);
        A(1:3:end,1:4) = xT1';
        A(1:3:end,13:16) = -(repmat(xT2(1,:),[4 1]).*xT1)';
        A(2:3:end,5:8) = xT1';
        A(2:3:end,13:16) = -(repmat(xT2(2,:),[4 1]).*xT1)';
        A(3:3:end,9:12) = xT1';
        A(3:3:end,13:16) = -(repmat(xT2(3,:),[4 1]).*xT1)';
      otherwise
        error('Number of dimensions not supported');
    end

    % (iii) (iv)
    [ disc disc V ] = svd( A );
    if size(x1,1)==2
      H = reshape(V(:,end),3,3)';
    else
      H = reshape(V(:,end),4,4)';
    end

    % Decondition
    H = (T2\H)*T1;
  case 10 % use RANSAC
    % Reference: HZ2, p.123, Algorithm 4.6
    nMax=0;
    for i=1:1000
      randSamp = randSample(size(x1,2),4);
      HTemp = computeHomography(x1(:,randSamp),x2(:,randSamp),1);
      x2Temp = normalizePoint( HTemp*normalizePoint( x1, -3 ), 3 );

      D = sum( ( x2Temp - x2  ).^2, 1 );

      nTemp = nnz( D<thres^2 );
      if nTemp>=nMax; nMax=nTemp; H=HTemp; end
    end
  case Inf % use PROSAC
end
