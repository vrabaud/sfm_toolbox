function demoCsfm( testCase )
% Simple demo file of CSFM on the shark data
%
% Just check the Rabaud Belongie CVPR09 paper
%
% The demos are as follows:
%  1: Shark data from Torresani et al, no noise
%  2: Orthographic absolute orientation computation
%  3: Orthographic exterior orientation computation
%  4: Orthographic rigid SFM examples
%
% USAGE
%  demoCSFM(testCase)
%  demoCSFM()
%
% INPUTS
%  testCase - number between 1 and 5
%
% OUTPUTS
%
% EXAMPLE
%
% See also 
%
% Vincent's Structure From Motion Toolbox      Version NEW
% Copyright (C) 2009 Vincent Rabaud.  [vrabaud-at-cs.ucsd.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the Lesser GPL [see external/lgpl.txt]

if(nargin<1), testCase = 1; end

% Generate data
switch testCase,
  case 1, % Torresani shark
    animGT=generateToyAnimation( 4.2, 'dt', 0 );
    animRank=5; animNBasis=2;
	animPeriod=120; % just for display purposes
  case 2, % artifical non-linear shape
    nFrame = 50; nPoint = 20;
    animGT=generateToyAnimation( 0.1,'nPoint',nPoint,'nFrame',nFrame,...
      'nBasis', 3, 'dt', 0, 'nLoop', 2 );
	animRank=9; animNBasis = 3;
	animPeriod=nFrame; % just for display purposes
  case 3, % artifical non-linear shape
    nFrame = 50; nPoint = 20;
    animGT=generateToyAnimation( 0.1,'nPoint',nPoint,'nFrame',nFrame,...
      'nBasis', [ 0 2 1], 'dt', 0, 'nLoop', 2 );
	animGT=animGT.addNoise('noiseS',0.01);
	animRank = 3; animNBasis = 3;
	animPeriod=nFrame; % just for display purposes
  case 4 % Walking Person
    animGT=generateToyAnimation( 5 );
	animRank = 6; animNBasis = 2;
	animPeriod=50; % just for display purposes
end

% pre-processing, precompute the optimal translation
anim=Animation(); anim.W=animGT.W;
anim.t=squeeze( mean(anim.W,2) ); anim.t(3,:)=0;

% Compute similarities within views
SimMin = computeAnimSimilarity(anim,-2);
imshow(SimMin)
disp('Displaying the similarities between frames');
disp('------------------------------------------------------------------');
in=input(['Press (n) to continue to next step, '...
  'anything else will quit the demo\n'],'s');
if in~='n'; return; end

% Compute samples to study with GNMDS
disp('Computing the samples to use in GNMDS (can take a while)');
minDist = 5;
samp = csfmComputeSample( anim, SimMin, 4000, minDist );
fprintf('\n');
disp('------------------------------------------------------------------');

% perform gnmds
disp('Computing GNMDS (can take a while too)');
[ lTot K ] = csfmGnmds( samp, SimMin, 240*8, 240*10, 1 );
anim.l=lTot(1:animNBasis,:);
disp('------------------------------------------------------------------');

% display the embedding for fun
if animNBasis==2 || animNBasis==3
	disp('Displaying some intermediate embedding maybe');
	switch animNBasis,
	case 2,
		plot(anim.l(1,:),anim.l(2,:),'.-');  hold on
	case 3,
		plot3(anim.l(1,:),anim.l(2,:),anim.l(3,:),'.-');  hold on
	end
	if animPeriod>0
	switch animNBasis,
		case 2,
		for i=1:anim.nFrame-animPeriod-1
			plot(anim.l(1,[i i+animPeriod]),anim.l(2,[i i+animPeriod]),'.-');
		end
		case 3,
		for i=1:anim.nFrame/2-1
			plot(anim.l(1,[i i+animPeriod]),anim.l(2,[i i+animPeriod]),'.-');
		end
	end
	end
	axis equal;
	fprintf('\n');
	disp(['--------------------------------------------------------', ...
	'----------']);
	in=input(['Press (n) to continue to next step, ', ...
	'anything else will quit the demo\n'],'s');
	if in~='n'; return; end
end

% Compute the best basis
% S is recomputed automatically
disp('Computing the resulting best SBasis');
SBasis=csfmComputeSBasis( anim.l, anim.W, animRank, 1, 1 );
save temp
anim.SBasis=SBasis;
% display the different elements of the basis
montageView(anim.SBasis);
disp('------------------------------------------------------------------');
in=input(['Press (n) to continue to next step, ' ...
'anything else will quit the demo\n'],'s');
if in~='n'; return; end

% Compute the best rotation matrices
anim.R = computeOrientation(anim.W,anim.S,'exteriorSequence');
disp('------------------------------------------------------------------');

% compute the error
[ err errFrame ]=anim.computeError('animGT',animGT);
fprintf( 'Average Reprojection Error: %f \n', mean(errFrame{1}(1,:)) );
fprintf( 'Average Depth Error: %f \n', mean(errFrame{1}(2,:)) );
fprintf( 'Average 3D Error: %f \n', err(1) );

% run the animation
playAnim(anim,'animGT',animGT,'showGT',true,'alignGT',true);

% perform gradient descent to optimize
disp('Perform the final optimization');
animIni=anim;
anim = bundleAdjustment( anim, 'nItrSBA', 64 );
disp('------------------------------------------------------------------');

% regenerate the 3D shapes
[ err errFrame ] = anim.computeError( 'animGT', animGT );
fprintf( 'Average Reprojection Error: %f \n', mean(errFrame{1}(1,:)) );
fprintf( 'Average Depth Error: %f \n', mean(errFrame{1}(2,:)) );
fprintf( 'Average 3D Error: %f \n', err(1) );

% display the final result close to the ground truth data
playAnim(anim,'animGT',animGT,'showGT',true,'alignGT',true);
