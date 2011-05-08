% Perform k-way normalized cut
%
% REFERENCE
%  Spectral Relaxation Models And Structure Analysis For K-Way Graph
%  Clustering And Bi-Clustering
%  Tech report, 2001
%  M. Gu, H. Zha, C. Ding, X. He, H. Simon, J. Xia
%
% USAGE
%  label = nCut( d, sigma, k )
%
% INPUTS
%  d           - distance matrix
%  sigma       -
%  k           - number of clusters
%
% OUTPUTS
%  label       -
%
% EXAMPLE
%
% See also
%
% Vincent's Structure From Motion Toolbox      Version 1.1
% Copyright (C) 2008-2011 Vincent Rabaud.  [vrabaud-at-cs.ucsd.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the GPL [see external/gpl.txt]

function label = nCut( d, sigma, k )

N = size( d, 1 );
W = exp(-d.^2/(2*sigma^2));

d = sum(W,1);

B = sparse(1:N,1:N,1./sqrt(d));
WHat = B*W*B; [ V D ] = eig(WHat);

% get the eigen vectors in descending order
[ disc dindex ] = sort(diag(D));
Yk = V(:,dindex(k:-1:1));

% find Xk
Yk = Yk./repmat(sqrt(sum(Yk.^2,2)),1,k);
label = kmeans2( Yk, k );
