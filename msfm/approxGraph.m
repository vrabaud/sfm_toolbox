function [ S2 S3Exp ] = approxGraph(S2,S3)
% Project a diadic/triadic graph to a 2D graph
%
%
% USAGE
%  [ S2 S3Exp ] = approxGraph(S2,S3)
%
% INPUTS
%  S2     - similarity matrixbetween pairs
%  S3     - similarity tensor between triplets
%
% OUTPUTS
%  S2     - new similarity graph
%  S3Exp  - S3 projected to 2D
%
% EXAMPLE
%
% See also GENERATETOYANIMATION, VIEWANIMSIMILARITY
%
% Vincent's Structure From Motion Toolbox      Version 3.0
% Copyright (C) 2008-2011 Vincent Rabaud.  [vrabaud-at-cs.ucsd.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the GPL [see external/gpl.txt]

ind=cell(nFrame,1);
for i=1:nFrame
  [j,k]=find(S3{i}(i+1:end,i+1:end));
  if isempty(j); continue; end
  temp=k<=j; j(temp)=[]; k(temp)=[];
  if isempty(j); continue; end
  ind{i}=[repmat(i,length(j),1), j+i,k+i];
end

% Replace each triadic connection by 3 diadic ones
M=cell2mat(ind);
S3Exp = zeros( nFrame, nFrame );
S3n = zeros( nFrame, nFrame );
for i=1:size(M,1)
  w=S3{M(i,1)}(M(i,2),M(i,3)); if w==0; w=Inf; end
  S3Exp(M(i,1),M(i,2)) = S3Exp(M(i,1),M(i,2)) + w;
  S3Exp(M(i,1),M(i,3)) = S3Exp(M(i,1),M(i,3)) + w;
  S3Exp(M(i,2),M(i,3)) = S3Exp(M(i,2),M(i,3)) + w;
  S3n(M(i,1),M(i,2)) = S3n(M(i,1),M(i,2)) + 1;
  S3n(M(i,1),M(i,3)) = S3n(M(i,1),M(i,3)) + 1;
  S3n(M(i,2),M(i,3)) = S3n(M(i,2),M(i,3)) + 1;
end

S3Exp = S3Exp + S3Exp';
S3n = S3n + S3n';
%S3Exp(S3Exp==0) = Inf;

%figure(1); imshow(-S2,[]);
%figure(2); imshow(-S3Exp,[]);

S2 = (1/2*S2 + 1/3*S3Exp)./( 1 + S3n );
for i=1:size(S2,1); S2(i,i) = 0; end

%figure(3); imshow(-S2,[]);
