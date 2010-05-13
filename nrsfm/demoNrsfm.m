function demoSFM( demoNumber )
% Demo for several abilities of LSML.
%
% Run with different integer values of demoNumber to run different demos.
%
% The demos are as follows:
%  1: Several visualizations
%  2: Orthographic NRSFM with random data
%  3: Orthographic NRSFM with random degenerate data
%  4: Orthographic NRSFM with shark data
%
% USAGE
%  demoSFM( demoNumber )
%
% INPUTS
%  demoNumber - [1] value between 1 and 4
%
% OUTPUTS
%
% EXAMPLE
%  demoSFM(1)
%
% See also
%
% Vincent's Structure From Motion Toolbox      Version 3.0
% Copyright (C) 2009 Vincent Rabaud.  [vrabaud-at-cs.ucsd.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the GPL [see external/gpl.txt]

if(nargin<1), demoNumber=1; end; c
if(demoNumber<1 || demoNumber>10), error('Invalid demo number.'); end
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
    if demoNumber<10; demoSFM(demoNumber+1); end
  case 'r'
    demoSFM(demoNumber);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function demo1() %#ok<DEFNU>
disp('Several Visualization Possibilities');

anim=generateToyAnimation(0.1,'nFrame',50,'nPoint',20,'nBasis',3,'nLoop',2);
disp('Playing an animation');
playAnim(anim);

disp('* Computing Similarities within the animation');
Sim = computeAnimSimilarity( anim, 2 );
disp('Visualize the similarities'); close all;
viewAnimSimilarity(anim,-Sim);
disp('Press enter to view ten random frames at the same time');
pause
close all;
montageView(anim.S(:,:,randSample(50,10)));
disp('Press enter for a person walking');
pause
close all;
anim=generateToyAnimation(5);
playAnim(anim, 'showConn', true, 'nCam', 100);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function demo2() %#ok<DEFNU>
disp('Orthographic NRSFM with random data');
% Test animations
nFrame = 50; nPoint = 20; nBasis = 3;
animGT=generateToyAnimation( 0.1,'nPoint',nPoint,'nFrame',nFrame,...
  'nBasis', nBasis, 'doSMean', 1, 'dt', 0 );
animGT=animGT.addNoise('noiseS',0);

%Xiao-Kanade
disp('**Computing NRSFM with Xiao''s method');
animXiao = computeNrsfm( 2, animGT.W );
[ err errFrame ] = animXiao.computeError( 'animGT', animGT );
fprintf( '\nXiao: Reprojection error %0.4f and 3D-error %0.4f\n\n', ...
  mean(errFrame{1}(1,:)),err(1));
playAnim(animXiao,'animGT',animGT,'showGT',true,'alignGT',true);

%Torresani
disp('**Computing NRSFM with Torresani''s method');
animTorr = computeNrsfm( 1, animGT.W, 'nBasis', 3, 'nItr', 50 );
[ err errFrame ] = animTorr.computeError( 'animGT', animGT );
fprintf( '\nTorresani: Reprojection error %0.4f and 3D-error %0.4f\n\n',...
  mean(errFrame{1}(1,:)),err(1));
playAnim(animTorr,'animGT',animGT,'showGT',true);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function demo3() %#ok<DEFNU>
disp('Orthographic NRSFM with random degenerate data');
%Xiao-Kanade with degenerate deformations
nFrame = 50; nPoint = 20; nBasis = [ 1 0 1 ];
animGT=generateToyAnimation( 0.1,'nPoint',nPoint,'nFrame',nFrame,...
  'nBasis', nBasis, 'dt', 0 );
animGT=animGT.addNoise('noiseS',0);
animXiao = computeNrsfm( 2, animGT.W );
[ err errFrame ] = animXiao.computeError( 'animGT', animGT );
fprintf( '\nXiao: Reprojection error %0.4f and 3D-error %0.4f\n\n', ...
  mean(errFrame{1}(1,:)),err(1));
playAnim(animXiao,'animGT',animGT,'showGT',true);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function demo4() %#ok<DEFNU>
disp('Orthographic NRSFM for the shark data');
%%% Shark data:
animGT=generateToyAnimation( 4.2,'dt', 0 );
animGT=animGT.addNoise('noiseS',0);
animXiao = computeNrsfm( 2, animGT.W);
[ errXiao errFrameXiao ] = animXiao.computeError( 'animGT', animGT );
fprintf( '\nXiao: Reprojection error %0.4f and 3D-error %0.4f\n\n', ...
  sum(errFrameXiao{1}(1,:)),errXiao(1));

animTorr = computeNrsfm( 1, animGT.W, 'nBasis', 3, 'nItr',30 );
[ errTorr errFrameTorr ] = animTorr.computeError( 'animGT', animGT );
fprintf( '\nTorresani: Reprojection error %0.4f and 3D-error %0.4f\n\n',...
  mean(errFrameTorr{1}(1,:)),errTorr(1));
%  playAnim(animXiao,'animGT',animGT,'showGT',true);


% perform gradient descent on the coefficients to improve the quality of
% the result
% animXiao = bundleAdjustment( animXiao, 'nItr', 10 );
% 
% [ errXiao errFrameXiao ] = animXiao.computeError( 'animGT', animGT );
% fprintf( 'Xiao+gradient descent: Reprojection error %0.4f and 3D-error %0.4f\n\n', ...
%   sum(errFrameXiao{1}(1,:)),errXiao(1));
% 
% playAnim(animXiao,'animGT',animGT,'showGT',true);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
