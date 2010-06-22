% Generate some shark data
animGT=generateToyAnimation( 4.2,'noise', 0.0, 'fixT', true );
dimSSpace = 2;
anim = computeNRSFM( 2, animGT.W);
[ err ] = computeNRSFMError( animGT, anim );
fprintf( 'Error with Xiao-Kanade: %e \n', sum(sum(err(1:2,:))) );
fprintf( 'Reprojection error: %e \n', computeSFMError( 'reproj', 0, 'anim', anim ) );

% Force some constraints
anim.SBasis(:,:,1)=-anim.SBasis(:,:,1);

anim.l(1,:)=1; anim.t(:)=0;
anim=anim.setFirstPRtToId();
anim = bundleAdjustment( anim, 'nItr', 100 );
anim=anim.setFirstPRtToId();

[ err ] = computeNRSFMError( animGT, anim );
fprintf( 'Error with Xiao-Kanade: %e \n', sum(sum(err(1:2,:))) );
fprintf( 'Reprojection error: %e \n', computeSFMError( 'reproj', 0, 'anim', anim ) );

% Set back to the original camera
for t=1:anim.nFrame
  anim.R(:,:,t)=animGT.R(:,:,t)'*anim.R(:,:,t);
end

animGT=generateToyAnimation( 4,'noise', 0.0, 'fixT', true );

anim.W=animGT.W;

[ err ] = computeNRSFMError( animGT, anim );
fprintf( 'Error with Xiao-Kanade: %e \n', sum(sum(err(1:2,:))) );
fprintf( 'Reprojection error: %e \n', computeSFMError( 'reproj', 0, 'anim', anim ) );

% Normalize the whole thingy
%  anim.SBasis(:,:,i)/mean(anim.l(1,:)); anim.l(1,:)=1;
for i=2:3
  anim.l(i,:)=norm(anim.SBasis(:,:,i),'fro')*anim.l(i,:);
  anim.SBasis(:,:,i)=anim.SBasis(:,:,i)/norm(anim.SBasis(:,:,i),'fro');
  anim.SBasis(:,:,1)=anim.SBasis(:,:,1)+mean(anim.l(i,:))*...
    anim.SBasis(:,:,i);
  anim.l(i,:)=anim.l(i,:)-mean(anim.l(i,:));
end

save('jawSource.mat','anim');