% Simple demo of MSFM
%
% MSFM is described in:
% Vincent Rabaud and Serge Belongie
% Re-Thinking Non-Rigid Structure From Motion
% CVPR 2008, Anchorage, Alaska.
%
% USAGE
%
% INPUTS
%
% OUTPUTS
%
% EXAMPLE
%
% See also GENERATETOYANIMATION, VIEWANIMSIMILARITY
%
% Vincent's Structure From Motion Toolbox      Version 3.0
% Copyright (C) 2009 Vincent Rabaud.  [vrabaud-at-cs.ucsd.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the GPL [see external/gpl.txt]

disp('------------------------------------------------------------------');
disp('Generating data');

nCase = 4.2;
switch nCase
  case 1
    %%% Bending Cylinder
    animGT = generateToyAnimation( 1, 'dt',1 );
  case 2
    %%% Roller Coaster
    animGT = generateToyAnimation( 2, 'dt',0 );
  case 3
    %%% Slinky
    animGT = generateToyAnimation( 3, 'dt', 1 );
  case 3.5
    %%% Real slinky
    animGT = generateToyAnimation(3.5);
  case 4.2
    %%% Shark
    animGT = generateToyAnimation( 4.2, 'dt', 0 );
end

% Compute pairwise similarities
disp('------------------------------------------------------------------');
disp('Computing pairwise similarities');
Sim2 = computeAnimSimilarity(animGT,2);
viewAnimSimilarity(animGT,-Sim2,'nCam',0);

% Compute triplet similarities
% Sim3 = computeAnimSimilarity(animGT,3,'Sim2',Sim2,'perc', 1);

% Compute the rigid shape chain
disp('------------------------------------------------------------------');
disp('Computing the rigid shape chain');
[animRSC,basis,ind] = msfmComputeRsc( animGT.W, Sim2);
%  [animRSC,basis,ind] = msfmComputeRsc( animGT, Sim2, Sim3);
playAnim(animRSC,'fps',20,'animGT',animGT);
cclf
montageView(anim.SBasis);

[ err errFrame errTot ] = animRSC.computeError( 'animGT', animGT );
fprintf( 'Average Reprojection Error: %f \n', mean(errFrame{1}(1,:)) );
fprintf( 'Average Depth Error: %f \n', mean(errFrame{1}(2,:)) );
fprintf( 'Average 3D Error: %f \n', err(1) );

% final optimization (actually, that willredo the two previous steps
% as they are very quick
disp('------------------------------------------------------------------');
disp('Computing the final optimizationrigid shape chain');
anim = nrsfmMsfm( animGT.W, 1, animGT );
