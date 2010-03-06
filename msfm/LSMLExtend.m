% X, XGt is Dxn
function XFinal = LSMLExtend( manifold, XGt, X, mask, xi, yi )

mask = logical( vect( mask, 'v' ) );

% get the closest ground truth point
temp = X(1:3:end,:);
lam = diff(minmax(temp(:)'))/20
done = repmat( false, [ 1 size(X,2) ]);
Xitr = inf( size( X ) );
for i=1:500
  sum(~done)
  if i==1
    [ dist ind ] = min( pdist2( X(mask,:)', XGt(mask,:)' ),[],2 );
    XitrDone = XGt;
  else
    XitrDone = Xitr(:,done);
    [ dist ind ] = min( pdist2( X(mask,~done)', XitrDone(mask,:)' ),[],2 );
  end
  
  toDo = dist<lam^2; temp = find(~done);
  ind = ind(toDo); toDo = temp(toDo);

  % Start from it and go to the sought point traversing the manifold
  H = LSMLcomputeH( XitrDone(:,ind), manifold, true );

  
%  i=20
%temp = LSMLcomputeH(XGt(3*20-2:3*20,:),manifold);
%temp = 30*temp(3*i-2:3*i,:,:);
%S = reshape(S,3,[],size(S,2));
%temp2 = squeeze(max(sum(temp.^2,1),[],2));
%temp(:,:,temp2>mean(temp2)+0.6*std(temp2))=0;
size(H)
size(XitrDone(3*20-2:3*20,ind))
visualizeManifold(squeeze(XitrDone(3*20-2:3*20,ind)),1,struct('VF',H(3*20-2:3*20,:,:)));


  pause
  
  
  
  % get the coeff that will make Xitr go from Xitr to X
  for j = 1:length(toDo)
    temp = X(mask,toDo(j))-XitrDone(mask,ind(j));
    tang = squeeze(H(:,:,j))*(H(mask,:,j)\temp);
    done(toDo(j))=true;
    if norm(tang)>10; 
      X(mask,toDo(j))
      XitrDone(mask,ind(j))
      H(mask,:,j)
      'coin', 
    continue
    end
    % should check for the length of the iteration 
    Xitr(:,toDo(j)) = XitrDone(:,ind(j)) + tang;
     
  end
%if rand()>0.999
  zi = reshape(Xitr(3*20,:),size(xi,1),size(xi,2));
  colormap gray; h=surf( xi, yi, zi, repmat( .7, size(zi) ) ); drawnow;
  shading flat; 
light('Position',[0 -2 1])
lighting phong; alpha( h, .6 );
    %end     
pause
  if nnz(~done)==0; break; end
end

XFinal = Xitr;
