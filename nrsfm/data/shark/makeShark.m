function A3=makeShark()

p = which('makeShark.m');
load( [ p(1:end-11) 'jaws' ] );

A3 = permute( reshape( P3_gt, [], 3, size(P3_gt,2) ), [2, 3, 1] );
