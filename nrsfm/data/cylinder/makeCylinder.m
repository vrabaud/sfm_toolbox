function [ A3 conn ] = makeCylinder( nPoint, nRound, nFrame )

% 0<=Phi<=pi, 0<=Th<=pi/2
Th=0; Phi=0; delTh=rand()/20; delPhi=rand()/5;

A3=zeros( 3, nPoint + (nPoint-2)*nRound, nFrame );
A3(:,1,:)=0;
rad=0.5/nPoint;
for i=1:nFrame
  % Compute the new positions of the points
  next=[ 0 0 1 ]'/nPoint;
  for j=2:nPoint
    next = rotationMatrix( cross( [ 0 0 1 ]',...
      [ cos(Phi) sin(Phi) 0 ]' ), Th*(1+2*j/nPoint) )*next;
    A3( :, j, i) = A3( :, j-1, i ) + next;
  end
  % create the points around
  for j=2:nPoint-1
    temp1 = cross( A3( :, j, i )-A3( :, j-1, i ), A3( :, j+1, i ) - ...
      A3( :, j, i ) );
    temp2 = A3( :, j-1, i )-A3( :, j, i ) + A3( :, j+1, i )-A3( :, j, i );
    
    if norm(temp1)~=0
      temp = -sin(Phi)*temp1/norm(temp1) + cos(Phi)*temp2/norm(temp2);
      temp = rad*temp/norm(temp);
    else
      temp1 = [ 1 0 0 ]; temp2 = [ 0 1 0 ]; temp=[ rad 0 0 ]';
    end
    
    for k=1:nRound
      A3( :, nPoint+nRound*(j-2)+k, i ) = A3( :, j, i ) + ...
        rotationMatrix( cross(temp1,temp2), 2*pi/nRound*(k-1) )*temp;
    end
  end
  
  Th = Th + delTh; Phi = Phi + delPhi;
  if Th<0 || Th>(pi/(nPoint+1)); delTh=-delTh; Th = Th + 2*delTh; end
  if Phi<0 || Phi>pi; delPhi=-delPhi; Phi = Phi + 2*delPhi; end
end

A3 = A3(:,nPoint+1:end,:);

% Create connectivity
conn = cell(1,1);
k=1;
for i=1:nPoint-2
  conn{k} = (i-1)*nRound + [ 1 : nRound, 1 ];
  k = k + 1;
end

% for i=2:nPoint-2
%   for j=1:nRound
%     k = k + 1; conn( k, : ) = [ j j+1 ];
%   end
%   conn( k, 2 ) = 1;
% end
