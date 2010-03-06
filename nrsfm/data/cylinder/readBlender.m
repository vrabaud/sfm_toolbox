addpath(genpath('~/matlab'));
addpath(genpath('./'));



anim=generateToyAnimation( struct('type',1) );

playAnimation( anim.A3, [], -20 )

playAnimation( anim.A2, [], -20 )
