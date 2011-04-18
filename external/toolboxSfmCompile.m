function toolboxSfmCompile(doSba,gccVer)
% Compiles all the private routines
%
% This function can recompile all that you need for the toolbox to run.
% That includes:
%  - the SBA static library that is needed to build the ...
%  - ... SBA mex file
%  - 3 normal mex files
% By default, the toolbox comes with precompiled binaries for windows 32
% It is easy to compile on Linux and I am developing on Linux so it works
% for sure there.
%
% USAGE
%  toolboxCompile()
%
% INPUTS
%  doSba      - [false] if true, it will compile the static library of SBA
%               as well as the mex file for SBA
%  gccVer     - on Linux, specify your gcc version here (e.g. 4.2 if
%               gcc-4.2 is a supported compiler on Matlab), which can be
%               useful if matlab thinks your gcc is too old/new
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

if nargin==0; doSba=false; end
if nargin<=1; gccVer=num2str(4.2); else gccVer=num2str(gccVer); end

if strcmp(computer,'GLNX86') || strcmp(computer,'GLNXA64')
  gcc_extra = [ 'CXX=g++-' gccVer ' CC=g++-' gccVer ' LD=g++-' gccVer ];
end

disp('Compiling.......................................');
savepwd=pwd; cd(fileparts(mfilename('fullpath'))); cd('../');

dir_external = [ pwd '/external/' ];
dir_sba = [ dir_external 'sba/'];

% build the SBA static library and the SBA mex file
if doSba
  cd(dir_sba);
  % delete previous object files
  fileList = {'sba_levmar.o', 'sba_levmar_wrap.o', 'sba_lapack.o', ...
    'sba_crsm.o', 'sba_chkjac.o', 'libsba.a', 'sba.lib', ...
    'matlab/sba.mexw32', 'matlab/sba.mexw64'};
  for i=1:length(fileList)
    if exist([dir_sba, '/', fileList{i}], 'file'), delete(fileList{i}); end
  end
  
  sba_files = [ 'sba.c ../sba_levmar_wrap.c ../sba_levmar.c ' ...
    '../sba_lapack.c ../sba_crsm.c ../sba_chkjac.c' ];
  
  cd matlab
  switch computer
    case {'PCWIN','PCWIN64'},
      ext = mexext;
      ext = ext(5:6);
      library_flags = [ '-D_WIN32 -I' dir_sba ' -L' dir_sba ...
        ' -L' dir_external 'blasLapack/win' ext ' -output sba.' mexext ' '];
      % Matlab on Windows 32/64
      if strcmp(computer, 'PCWIN')
        eval([ 'mex ' library_flags ...
          '-lclapack.lib -lblas.lib -llibF77.lib ' ...
          '-llibI77.lib' sba_files ]);
      else
        eval([ 'mex LINKFLAGS="$LINKFLAGS /NODEFAULTLIB:libcmt.lib" ' ...
          library_flags '-lclapack_nowrap -lf2c -lBLAS_nowrap ' sba_files ]);
      end
    case {'GLNX86'},
      % Matlab on Linux
      eval([ 'mex ' gcc_extra ' -I../ -O -ldl ' sba_files ...
        '/usr/lib/liblapack.a /usr/lib/atlas/libblas.a ' ...
        '/usr/lib/libgfortran.so.3' ]);
    case {'GLNXA64'},
      % Matlab on Linux 64
      eval([ 'mex -lgfortran ' gcc_extra ' -I../ -O -ldl ' sba_files ...
        ' /usr/lib/liblapack_pic.a ' ...
        '../../blasLapack/linux64/libf77blas.a ' ...
        '../../blasLapack/linux64/libatlas.a '
        ]);
    case {'i686-pc-linux-gnu', 'i686-unknown-linux-gnu', ...
        'x86_64-pc-linux-gnu', 'x86_64-unknown-linux-gnu'},
      % Octave on Linux
      mkoctfile --mex ./sba.c -I../ -lsba -L../
  end
  
  cd ../../..
end

cd sfm/private/sba
switch computer
  case {'PCWIN','PCWIN64'},
    if doSba
      % Matlab on Windows 32/64
      system('cl /nologo /O2 sbaProjection.c /link /dll /out:sbaProjection.dll');
    end
  case {'GLNX86','GLNXA64','i686-pc-linux-gnu', 'i686-unknown-linux-gnu', ...
      'x86_64-pc-linux-gnu', 'x86_64-unknown-linux-gnu'},
    % Matlab and Octave on Linux
    system('gcc -Wall -fPIC -O3 -shared -o sbaProjection.so sbaProjection.c');
end
cd ../../..

% compile the remaining mex files
rd=fileparts(mfilename('fullpath')); rd=rd(1:end-9);

% general compile options (can make architecture specific)
optsAfter={};
switch computer
  case {'PCWIN','PCWIN64'},
    if strcmp(computer, 'PCWIN')
      win_folder = 'win32';
    else
      win_folder = 'win64';
    end
    % Matlab on Windows
    lapacklib = fullfile(matlabroot, 'extern', 'lib', win_folder, ...
      'microsoft', 'libmwlapack.lib');
    blaslib = fullfile(matlabroot, 'extern', 'lib', win_folder, ...
      'microsoft', 'libmwblas.lib');
    opts={'-output'};
    optsAfter = {lapacklib, blaslib};
  case {'GLNX86', 'GLNXA64'},
    % Matlab on Linux
    % if you get warnings on linux, you could force the gcc version by
    % adding those options: 'CXX=g++-4.1' 'CC=g++-4.1' 'LD=g++-4.1'
    opts = { ['CXX=g++-' gccVer] ['CC=g++-' gccVer] ['LD=g++-' ...
      gccVer] '-largeArrayDims' '-l' 'mwlapack' '-l' 'mwblas' '-output'};
  case {'i686-pc-linux-gnu', 'i686-unknown-linux-gnu', ...
      'x86_64-pc-linux-gnu', 'x86_64-unknown-linux-gnu'},
    % Octave on Linux
    opts = {'-o'};
end

% general compile options
fs={'computeCsfmInfimumMex','refineExteriorOrientationMex'};%, ...
%  'msfmRotationDerivativeMex'};
ds={'nrsfm', 'sfm', 'msfm'};
for i=1:length(fs)
  mex(opts{:},[rd '/' ds{i} '/private/' fs{i} '.' mexext], ...
    [rd '/' ds{i} '/private/' fs{i} '.c'], optsAfter{:});
end

cd(savepwd); disp('..................................Done Compiling');

end
