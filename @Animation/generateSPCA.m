function anim = generateSPCA( anim, nPCA )
% Generate the 3D PCA approximation of an Animation
%
% It simply performs PCA on S
%
% USAGE
%  anim = anim.generateSPCA()
%
% INPUTS
%  anim     - Animation object (help Animation for details)
%
% OUTPUTS
%  anim     - modified Animation
%
% EXAMPLE
%
% See also GENERATETOYANIMATION
%
% Vincent's Structure From Motion Toolbox      Version NEW
% Copyright (C) 2009 Vincent Rabaud.  [vrabaud-at-cs.ucsd.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the Lesser GPL [see external/lgpl.txt]

if isempty(anim.S); error('No S in the Animation'); end

[ U mu vars ]=pca(anim.S(:,:,:));

anim.l=[]; anim.SBasis=[];
anim.l = [ ones( 1, anim.nFrame ); Yk ];
anim.SBasis=reshape(U(:,1:nPCA),3,anim.nPoint,nPCA);
anim.SBasis(:,:,end+1)=mu;
anim.SBasis=anim.SBasis(:,:, [ end, 1:end-1 ] );
