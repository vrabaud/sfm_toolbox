function S = makeSlinky(len,radius,nFrame,caseN)

if nargin==4 % it is the real data case
  load( [ '/daVinci/slinky/mat' int2str(caseN) '/res01500.mat' ] );
  k=1; traj = zeros(3,2,3);
  for j=1:length(cen{1}) %#ok<USENS>
    cenPrev = cen{1}(j,:);
    if caseN==1 && cenPrev(2)>410; continue; end %1
    if caseN==2 && cenPrev(2)>400 || (cenPrev(2)>310 && cenPrev(1)<200); continue; end %2
    if caseN==3 && cenPrev(2)>260 && cenPrev(1)<340; continue; end %3
    if caseN==4 && cenPrev(2)>450; continue; end %3
    traj(1,:,k) = cenPrev;
    for i=2:1500
      
      %plot(cen{i}(:,1),cen{i}(:,2),'.'); drawnow;
      %end
      
      [ disc ind ] = min( pdist2( traj(i-1,:,k), cen{i}) );
      traj(i,:,k) = cen{i}(ind,:);
    end
    k=k+1;
  end
  for i=1:1500
    temp = pdist( squeeze(traj(i,:,:))' );
    if min( temp(:) )<=eps
      temp = squareform(temp);
      temp = temp + diag(Inf*ones(1,size(temp,2)));
      [y x] = find(temp<=eps);
      x = unique( [ x' y' ]);

      temp = pdist2( traj(i,:,x(1)), squeeze(traj(i-1,:,:))' );
      [ disc ind ] = min(temp);
      traj(:,:,setdiff( x, ind )) = [];
    end
  end

  %if caseN==2; traj = traj(100:end,:,:); end
%   temp=anim.W(:,:,100)-anim.W(:,:,101);
%   temp=anim.W(:,:,1:end-1)-anim.W(:,:,2:end);
%   temp=sqrt(temp(1,:).^2+temp(2,:).^2);
%   plot(temp)
  
  
  anim.W = permute( traj, [2,3,1]);
  anim.W(2,:) = -anim.W(2,:);
  anim.W = anim.W(:,:,1:5:end);
  anim.isProj = false;
  S = anim;
else
  nTraj=10;
  for i=1:nTraj
    rad=1+radius*2*rand();
    deltaT=len/3000; deltath=pi/500; th=-pi/2-deltath; j=1;
    A{i}=zeros(3,100);
    for T=-len/2:deltaT:len/2; A{i}(:,j)=[T rad 1]'; j=j+1; end
    th=pi/2-deltath;
    while th>-pi/2
      A{i}(:,j)=[len/2+rad*cos(th) rad*sin(th) 1]'; th=th-deltath; j=j+1;
    end
    for T=len/2:-deltaT:-len/2; A{i}(:,j)=[T -rad 1]'; j=j+1; end
    A{i}(3,:)=sqrt( radius^2 - (rad-1-radius)^2);

    A{i}=A{i}(:,1:j-1);
  end


  % Make it evolve through time
  nPoint=size(A{1},2);
  speed=20; nSamp=100; first=1;
  nPointSamp=length(1:nSamp:nPoint)-24;
  sense=1;

  for i=1:nFrame
    inter=first:nSamp:first+nPointSamp*nSamp;

    if inter(end)>size(A{1},2) || inter(1)<1
      sense=-sense;
    end
    while inter(end)>size(A{1},2) || inter(1)<1
      first=first+sense*speed;
      inter=first:nSamp:first+nPointSamp*nSamp;
    end

    B=zeros(3,0);
    for j=1:nTraj
      B=[B,A{j}(:,inter)];
    end

    first=first+sense*speed;

    if i==1; S=repmat(B,[1,1,nFrame]); continue; end
    S(:,:,i)=B;
    if rand()<0.1; speed=round(10+rand()*100); end
  end

  S = S( [2 3 1], :, : );
end
