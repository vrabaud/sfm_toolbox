% Compute manifold warping function.
%
% USAGE
%   manifold = LSMLoptimizeTH( X, Nmat, pw )
%
% INPUT
%   X or X3   - data [D x n] or [nSet x nPoint x D] or matrix
%   Nmat      - [n x n] neighborhood matrix [see computeNeighbor]
%   pw        - parameters struct (all but d are optional)
%     [parameters for rbfs]
%       .rbfFlag    - [1]   flag for using rbfs
%       .rbfscale   - [2]   [see rbfComputeBasis]
%       .rbfcluster - [1]   [see rbfComputeBasis]
%       .rbfk       - []    [see rbfComputeBasis]
%     [parameters for coordinate descent]
%       .d          - [REQ]  manifold dimension
%       .nrestart   - [25]   number of restarts
%       .nitr       - [50]   number of iterations
%       .lam        - [1e-3] regularization term
%       .THinit     - [rand] initial value for TH
%       .nSamples   - [1000] num to subsample for training
%       .show       - [0]    optionally display progress in figure(show)
%       .percTest   - [20]   percentage of testing data to consider
%
% OUTPUT
%   manifold
%     [dimensions]
%       .D          - original dimension
%       .f          - number of features per point
%       .d          - manifold dimension
%     [rbf]
%       .rbfFlag    - 1 if radial basis functions are being used
%       .rbfBasis   - basis functions ([] if rbfFlag==0)
%     [learned warp]
%       .pw         - parameters used to compute warping
%       .TH         - [d x f x D] computed warping function
%       .THerr      - error of computed warping function
%       .THerrF     - non-regularized error of computed warping function
%       .errTest    - test error

function manifoldBest = fwLSMLoptimizeTH( X, Nmat, pw )

nItr=pw.nrestart; pw.nrestart=1;

% system('fusermount -u /projects/fwHome');
% system('fusermount -u /projects/fwScratch');

system('sshfs vrabaud@137.110.131.21:/home/vrabaud /projects/fwHome');
system(['sshfs vrabaud@137.110.131.21:' ...
  '/scratch/vrabaud/ /projects/fwScratch']);

while 1
  try
    save('/projects/fwHome/SSFM/data','X','Nmat','pw');
  catch
    pause(1); continue;
  end
  break;
end

scriptName = fwLaunch( nItr, 'fwLSMLoptimizeTHNode', '~/SSFM/', ...
  '~/tempRes/', '~/temp/' );

nDone=0; pathToCheck='/projects/fwHome/tempRes/';
errBest=Inf;
while nDone<nItr-5
  name = dir(pathToCheck);
  pause(5); if isempty(name); continue; end
  for i = 1 : length( name )
    if isempty( strfind( name(i).name, scriptName ) ); continue; end
    fileName=[ pathToCheck '/' name(i).name ];
    %if isempty(regexp( fileName, fileId )); continue; end
    %if any(strcmp(fileName, fileDone)); continue; end
      
    % In case the file is being written
    try load( fileName ); catch continue; end
    if manifold.THerr<errBest
      errBest = manifold.THerr;
      manifoldBest = manifold;
    end
    system( sprintf( 'rm %s', fileName ) );
    %system( sprintf( 'rm %s', ['/projects/fwHome/SSFM/' name(i).name '*' ]) );
    nDone = nDone + 1;
    %nDone
  end
end

%pw.THinit = manifoldBest.TH;
%manifold = LSMLoptimizeTH( X, Nmat, pw );
