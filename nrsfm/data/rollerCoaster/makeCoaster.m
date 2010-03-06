function S=makeCoaster(len,rad,nFrame,doSymmetric)

deltaT=len/1000; deltath=pi/1000; th=-pi/2-deltath; j=1; A=zeros(3,100);
while th>-3*pi/2
  A(:,j)=[-len/2+rad*cos(th) rad*sin(th) 1]'; th=th-deltath; j=j+1;
end
for T=-len/2:deltaT:len/2; A(:,j)=[T rad 1]'; j=j+1; end
th=pi/2-deltath;
while th>-pi/2
  A(:,j)=[len/2+rad*cos(th) rad*sin(th) 1]'; th=th-deltath; j=j+1;
end
for T=len/2:-deltaT:-len/2; A(:,j)=[T -rad 1]'; j=j+1; end
A=A(:,1:j-1);

% Bend the bike chain
nPoint=size(A,2);
if doSymmetric
  center=0;
else
  center=len/3;
end
A(3,:)=6*(0.5-abs((A(1,:)-center)/(len/2+rad)));

% Make it evolve over time
speed=50; nSamp=100; first=1;
nPointSamp=length(1:nSamp:nPoint);

for i=1:nFrame
  inter=first:nSamp:first+nPointSamp*nSamp;
  inter=mod(inter+speed-1,nPoint)+1;
  B=A(:,inter);
  first=first+speed;
  if i==1; S=repmat(B,[1,1,nFrame]); continue; end
  S(:,:,i)=B;
  if rand()<0.1; speed=round(50+rand()*50); end
end
