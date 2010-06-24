function demoSfm( demoNumber )
% Demo for several abilities of LSML.
%
% Run with different integer values of demoNumber to run different demos.
%
% The demos are as follows:
%  1: Homography computation with or without outliers
%  2: Homography computation with bundle adjustment
%  3: Orthographic absolute orientation computation
%  4: Orthographic exterior orientation computation
%  5: Orthographic rigid SFM examples
%  6: Bundle adjustment on camera parameters
%  7: Projective rigid SFM examples
%
% USAGE
%  demoSfm( demoNumber )
%
% INPUTS
%  demoNumber - [1] value between 1 and 7
%
% OUTPUTS
%
% EXAMPLE
%  demoSfm(1)
%
% See also
%
% Vincent's Structure From Motion Toolbox      Version 3.0
% Copyright (C) 2009 Vincent Rabaud.  [vrabaud-at-cs.ucsd.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the GPL [see external/gpl.txt]

if(nargin<1), demoNumber=1; end; c
if(demoNumber<1 || demoNumber>7), error('Invalid demo number.'); end
disp('Demos of various functions in SFM toolbox.');
disp('Steps marked w ''*'' or ''**'' may take time, please be patient.');
disp('------------------------------------------------------------------');
fprintf('BEGIN DEMO %i: ',demoNumber);
eval(['demo' int2str(demoNumber) '()']);
disp(['END OF DEMO ' int2str(demoNumber) '.']);
disp('------------------------------------------------------------------');
in=input(['Press (n) to continue to next demo or (r) to repeat demo\n'...
  'Anything else will quit the demo\n'],'s');
switch in
  case 'n'
    if demoNumber<7; demoSfm(demoNumber+1); end
  case 'r'
    demoSfm(demoNumber);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function demo1() %#ok<DEFNU>
disp('Homography example');
disp('Create 100 data points and a random homography');
x1 = rand( 2, 100 ); H0 = rand( 3, 3 );
x2 = normalizePoint( H0*normalizePoint( x1, -3 ), 3 );
H1 = computeHomography( x1, x2, 1 );
disp('Ratio of recovered to original homography');
H1./H0 %#ok<NOPRT>

% Inclusion of outliers: use RANSAC
disp('Creating 20%% of outliers');
x2(:,81:100) = rand( 2, 20 );
[ H2 n ] = computeHomography( x1, x2, 10, 0.01 );

fprintf( 'Recovered a homography with %d%% inliers\n', n );
disp('Ratio of recovered to original homography');
H2./H0 %#ok<NOPRT>
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function demo2() %#ok<DEFNU>
disp('Homography with bundle adjustment example');
disp('Create points and random view');

anim=Animation(); anim.W=zeros(2,100,5);
X=rand( 2, 100 );
for i=1:5
  H = rand( 3, 3 ); anim.mask(:,i)=logical(rand(100,1)>0.5);
  if i>=1; H=eye(3); end
  anim.W(:,anim.mask(:,i),i) = normalizePoint( H*normalizePoint( X(:,anim.mask(:,i)), -3 ), 3 );
end

% add some noise to the measurements
anim=anim.addNoise('noiseW', 0.01, 'doFillW', true);

% compute best homographies between those points
P=zeros(3,3,anim.nFrame);
anim.S=X;
for i=1:5
  mask=anim.mask(:,1) & anim.mask(:,i);
  P(:,:,i)=computeHomography( anim.W(:,mask,1), anim.W(:,mask,i), 1 );
end
anim.P=P;

err = anim.computeError();
errTot = [ err(1), 0 ];

[ anim ] = bundleAdjustment( anim );

err = anim.computeError();
errTot(2) = err(1);

out =sprintf( 'Reprojection error %0.4f/%0.4f, before/after BA\n\n',...
  errTot(1), errTot(2) );
fprintf(out);
% color=colorcube(8);
color = [ 1 1 0; 0 1 1; 1 0 1; 1 0 0; 0 1 0; 0 0 1; 0 0 0; 1 1 1 ];

for i=1:5
  plot(anim.W(1,anim.mask(:,i),i),anim.W(2,anim.mask(:,i),i),'.','color',...
    color(i,:));
  hold on;
end
plot(anim.S(1,:),anim.S(2,:),'ro'); axis tight;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function demo3() %#ok<DEFNU>
disp('Orthographic absolute orientation computation');
disp('Generate some random rigid data');
% Test animations
nFrame = 10; nPoint = 50;
animGT=generateToyAnimation( 0,'nPoint',nPoint,'nFrame',nFrame);

S1 = animGT.S+randn(3,animGT.nPoint)*0.05;
S2 = animGT.S+randn(3,animGT.nPoint)*0.05;

tt=randSample(animGT.nFrame,2);
[R t s] = computeOrientation(animGT.R(:,:,tt(1))*S1+...
  animGT.t(:,tt(1)*ones(1,animGT.nPoint)),S2, 'absolute');

disp([ 'The absolute orientation between two noisy 3D shapes gives a ' ...
  'difference of:']);
disp( 'Rotation Difference');
R-animGT.R(:,:,tt(1)) %#ok<NOPRT,MNEFF>
disp( 'Translation Difference');
t-animGT.t(:,tt(1)) %#ok<NOPRT,MNEFF>
disp( 'Scale Difference');
s-1 %#ok<NOPRT,MNEFF>
disp( '(Overall close to identity)');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function demo4() %#ok<DEFNU>
disp('Orthographic exterior orientation computation');
disp('Generate some random rigid data');
% Test animations
nFrame = 10; nPoint = 50;
animGT=generateToyAnimation( 0,'nPoint',nPoint,'nFrame',nFrame);
animGT=animGT.addNoise('noiseS', '5', 'doFillW', true);
disp('*Compute the exterior orientation using Gloptipoly 3');
t0=randSample(animGT.nFrame,1);

[R t] = computeOrientation( animGT.S(:,:,1),animGT.W(:,:,t0),...
  'exterior');

disp(['The exterior orientation between a shape reconstructed from '...
  'noisy data and a reprojected frame gives a difference of:']);
disp( 'Rotation Difference');
R-animGT.R(:,:,t0) %#ok<NOPRT,MNEFF>
disp( 'Translation Difference');
t(1:2)-animGT.t(1:2,t0) %#ok<NOPRT,MNEFF>
%  disp( 'Scale');
%  s-1
disp( '(Overall close to identity)');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function demo5() %#ok<DEFNU>
disp('Orthographic rigid SFM examples');
disp('Generate some random rigid data');
% Test animations
nFrame = 10; nPoint = 50;
animGT=generateToyAnimation( 0,'nPoint',nPoint,'nFrame',nFrame);
animGT=animGT.addNoise('noiseS','5','doFillW',true);

if ~exist('OCTAVE_VERSION','builtin')
  playAnim( animGT, 'frame', 1, 'nCam', 20 );
end

% The following should give the same low errors
fprintf( [ 'Computing several reconstructions on %d frames with %d '...
  'noisy features.\n' ], nFrame, nPoint );
anim = cell(5,2);
err = zeros(5,2); err3D = err;

for i = 1 : 5
  for j = 1 : 2
    switch i
      case 1,
        animGTSample=animGT.sampleFrame(1:3);
        anim{1,j} = computeSMFromW( false, animGTSample.W, 'method', 0, ...
          'nItrSBA', (j-1)*100 );
        out = '3 views, uncalibrated cameras:\n';
        typeTransform = 'homography';
        type3DError = 'up to a projective transform';
      case 2,
        animGTSample=animGT;
        anim{2,j} = computeSMFromW( false, animGT.W, 'method', 0, ...
          'nItrSBA', (j-1)*100 );
        out = 'All the views, uncalibrated cameras:\n';
        typeTransform = 'homography';
        type3DError = 'up to a projective transform';
      case 3,
        animGTSample=animGT.sampleFrame(1:3);
        anim{3,j} = computeSMFromW( false, animGTSample.W, 'method', 0, ...
          'isCalibrated', true, 'nItrSBA', (j-1)*100 );
        out = '3 views, calibrated cameras:\n';
        typeTransform = 'camera';
        type3DError = '';
      case 4,
        animGTSample=animGT;
        anim{4,j} = computeSMFromW( false, animGT.W, 'method', 0, ...
          'isCalibrated', true, 'nItrSBA', (j-1)*100 );
        out = 'All the views, calibrated cameras:\n';
        typeTransform = 'camera';
        type3DError = '';
      case 5,
        animGTSample=animGT;
        anim{5,j} = computeSMFromW( false, animGT.W, 'method', inf, ...
          'isCalibrated', true, 'nItrSBA', (j-1)*100 );
        out = 'All the views, calibrated cameras with SDP solving:\n';
        typeTransform = 'camera';
        type3DError = '';
    end
    errTmp = anim{i,j}.computeError();
    err(i,j) = errTmp(1);
    err3DTmp = anim{i,j}.computeError('animGT', animGTSample, ...
      'checkTransform', typeTransform );
    err3D(i,j) = err3DTmp(1);
  end
  out =sprintf([ out 'Reprojection error %0.4f/%0.4f and 3D error %s %0.4f/%0.4f , before/after BA\n\n'],...
    err(i,1), err(i,2), type3DError, err3D(i,1), err3D(i,2) );
  fprintf(out);
end

disp('The following keys can be used during an animation');
disp('Type ''help playAnim'' to display the following');
disp('''c'' toggle the connectivity display');
disp('''f'' toggle the first point display');
disp('''g'' toggle the ground truth display');
disp('''p'' set the axes to be pretty');
disp('''q'' quit the playing of the animation');
disp('''r'' rotate the axes');
disp('''t'' toggle the title display');
disp( [ '{''1'',''2'',''3''} switch the camera view ' ...
  '(1: fixed axes, 2: camera projection, 3:fixed camera}']);
disp('''numpad0'' set the view back to the original');
disp('''-'' slow down the speed');
disp('''+'' increase the speed');
disp('''space'' pause the animation');
disp([ '{''leftarrow'', ''rightarrow'', ''uparrow'', ''downarrow''} ' ...
  'rotate the 3D view']);

% Plot in 3D
animAligned=anim{4,2};
[ disc, disc, disc, disc, animAligned.S ] = animAligned.alignTo( ...
  animGT, 'camera');
animAligned.R=repmat(eye(3),[1,1,animAligned.nFrame]);
animAligned.t=zeros(3,animAligned.nFrame);
playAnim(animAligned,'animGT',animGT,'nCam',-1,'showGT',true,'alignGT',true);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function demo6() %#ok<DEFNU>
disp('Simple bundle adjustment with calibration parameters');
disp('Generate some random rigid data');
% Test animations
nFrame = 5; nPoint = 50; isProj = true;
animGT=generateToyAnimation( 0,'nPoint',nPoint,'nFrame',nFrame,...
  'isProj', isProj );
animGT.K = 5*(rand(3,nFrame)); animGT.K(2,:)=0;
animGT.W=animGT.generateW();

% Add some noise to this ground truth animation
S0 = animGT.S + randn(3,nPoint)*0.1/max(abs(animGT.S(:)));
R0 = quaternion(quaternion(animGT.R)+randn(4,nFrame)*0.001);
t0 = animGT.t+randn(3,nFrame)*0.1*max(abs(animGT.t(:)));
K0 = animGT.K+randn(3,nFrame)*0.001*max(abs(animGT.K(:))); K0(2,:)=0;

anim=Animation; anim.R=R0; anim.t=t0; anim.K=K0; anim.S=S0; anim.W=animGT.W;
errTmp = anim.computeError();
err3DTmp = anim.computeError('animGT', animGT, ...
  'checkTransform', 'rigid' );

err(1) = errTmp(1); err3D(1) = err3DTmp(1);

% perform bundle adjutment
[ anim info ] = bundleAdjustment( anim, ...
  'KMask', [ 0 1 0 ], 'nFrameFixed', 0 );
anim=anim.setFirstPRtToId();

errTmp = anim.computeError();
err3DTmp = anim.computeError('animGT', animGT, 'checkTransform', 'rigid' );
err(2) = errTmp(1); err3D(2) = err3DTmp(1);

out =sprintf([ 'Reprojection error %0.4f/%0.4f and 3D error %0.4f/%0.4f'...
  ', before/after BA\n\n'], err(1), err(2), err3D(1), err3D(2) );
fprintf(out);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function demo7() %#ok<DEFNU>
%%% Check for projective camera
disp('Projective rigid SFM examples');
% Test animations
nFrame = 10; nPoint = 50;
animGT=generateToyAnimation( 0,'nPoint',nPoint,'nFrame',nFrame,...
  'isProj',true,'dR', 1 );
animGT=animGT.addNoise('noiseS', '5', 'doFillW', true);
playAnim( animGT, 'frame', 1, 'nCam', 20 );
save('badAffine');
% The following should give the same low errors
fprintf( ['Computing several reconstructions on %d frames with %d noisy'...
  ' features.\n' ], nFrame, nPoint );
anim = cell(5,2);
err = zeros(5,2); err3D = err;

for i = 1 : 5
  for j = 1 : 2
    switch i
      case 1,
        animGTSample=animGT.sampleFrame(1:2);
        anim{1,j} = computeSMFromW( true, ...
          animGTSample.W, 'method', 0,'nItrSBA', (j-1)*100 );
        out = '2 views, uncalibrated cameras:\n';
        typeTransform = 'homography';
        type3DError = 'up to a projective transform';
      case 2,
        animGTSample=animGT.sampleFrame(1:2);
        anim{2,j} = computeSMFromW( true, ...
          animGTSample.W, 'method', 0,'isCalibrated',true,...
          'nItrSBA', (j-1)*100);
        out = '2 views, calibrated cameras:\n';
        typeTransform = 'homography';
        type3DError = '';
      case 3,
        animGTSample=animGT;
        anim{3,j} = computeSMFromW( true, ...
          animGT.W, 'method', 0, 'nItrSBA', (j-1)*100 );
        out = 'All the views, uncalibrated cameras, Sturm Triggs:\n';
        typeTransform = 'homography';
        type3DError = 'up to a projective transform';
      case 4,
        animGTSample=animGT;
        anim{4,j} = computeSMFromW( true, ...
          animGT.W, 'method', Inf, 'nItrSBA', (j-1)*100 );
        out = 'All the views, uncalibrated cameras, Oliensis Hartley:\n';
        typeTransform = 'homography';
        type3DError = 'up to a projective transform';
      case 5,
        animGTSample=animGT;
        anim{5,j} = computeSMFromW( true, ...
          animGT.W, 'method', Inf, 'isCalibrated',true,...
          'nItrSBA', (j-1)*100, 'doAffineUpgrade', true, 'nItrAff', 20 );
        typeTransform = 'rigid+scale';
        type3DError = 'up to a scaled rigid transform';
        out = 'All the views, calibrated cameras:\n';
    end
    errTmp = anim{i,j}.computeError();
    err(i,j) = errTmp(1);
    err3DTmp = anim{i,j}.computeError('animGT', animGTSample, ...
      'checkTransform',typeTransform);
    err3D(i,j) = err3DTmp(1);
  end
  out =sprintf([ out 'Reprojection error %0.4f/%0.4f and 3D-error %s %0.4f/%0.4f , before/after BA\n\n'],...
    err(i,1), err(i,2), type3DError, err3D(i,1), err3D(i,2) );
  fprintf(out);
end

% Plot in 3D
%  [ disc disc disc SBest ] = computeSFMError('3D', false, 'SGT', animGT.S(:,:,1), ...
%       'S', S(:,:,3,2), 'checkTransform', 'rigid' );
%
%  playAnim(anim,'animGT',animGT,'nCam',-1,'showGT',true,'alignGT',true);

end
