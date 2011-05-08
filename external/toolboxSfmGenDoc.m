function toolboxSfmGenDoc()
% Generate documentation, must run from dir toolbox.
%
% Requires external/m2html to be in path.
%
% Prior to running, update a few things in the overview.html file in:
%   toolbox/external/m2html/templates/frame-vincent/overview.html
%   1) The version / date
%   2) Link to rar/zip file
% Might be useful to use:
%
% USAGE
%  toolboxGenDoc
%
% INPUTS
%
% OUTPUTS
%
% EXAMPLE
%
% See also
%
% Vincent's Structure From Motion Toolbox      Version 3.0
% Copyright (C) 2008-2011 Vincent Rabaud.  [vrabaud-at-cs.ucsd.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the GPL [see external/gpl.txt]

% run m2html
params = {'mfiles',{'@Animation','csfm','linearAlgebra','nrsfm','sfm','visualize'}};
params = {params{:},'htmldir','doc','recursive','on','source','off'};
matlabpathTmp = matlabpath;
tmp = which('m2html','-all');
for i = 1 : length(tmp); rmpath( fileparts( tmp{i} ) ); end
addpath( genpath( fileparts( mfilename('fullpath') ) ) );
params = {params{:},'template','frame-vincent','index','menu','global','on'};
m2html(params{:});
addpath( matlabpathTmp );

% copy custom menu.html
copyfile('external\m2html\templates\menu-for-frame-vincent.html','doc/menu.html');

% copy history file
copyfile('external\history.txt','doc/history.txt');

% remove links to .svn and private in the menu.html files
d = { './doc' };

while ~isempty(d)
  dTmp = dir(d{1});
  for i = 1 : length(dTmp)
    name = dTmp(i).name;
    if strcmp( name,'.') || strcmp( name,'..'); continue; end
    if dTmp(i).isdir; d{end+1} = [ d{1} '/' name ]; continue; end
    if ~strcmp( name,'menu.html'); continue; end
    fid = fopen( [ d{1} '/' name ], 'r' ); c = fread(fid, '*char')'; fclose( fid );
    c = regexprep( c, '<li>([^<]*[<]?[^<]*)\.svn([^<]*[<]?[^<]*)</li>', '');
    c = regexprep( c, '<li>([^<]*[<]?[^<]*)private([^<]*[<]?[^<]*)</li>', '');
    fid = fopen( [ d{1} '/' name ], 'w' ); fwrite( fid, c ); fclose( fid );
  end
  d(1) = [];
end

% Zip the whole toolbox
s = dir( './' );
goodFolder = { '@Animation' 'csfm' 'doc' 'external' 'linearAlgebra' 'nrsfm' 'sfm' 'visualize' };
for i = length(s) : -1 : 1
  isGood = false;
  for j = 1 : length(goodFolder)
    if strcmp( s(i).name, goodFolder{j} )
      isGood = true;
    end
  end
  if isGood
    s(i).name = [ './' s(i).name ];
  else s(i) = [];
  end
end
files = cell(1,0);

while ~isempty(s)
  if s(1).isdir
    ls = length(s);
    s = [ s; dir( s(1).name ) ];
    for i = length(s) : -1 : ls+1
      if ~strcmp( s(i).name(1), '.' )
        s(i).name = [ s(1).name '/' s(i).name ];
      else s(i) = [];
      end
    end
    s(1) = [];
  else
    if ~strcmp( s(1).name(end), '~' )
      files{end+1} = s(1).name;
    end
    s(1) = [];
  end
end

zip( 'vincentToolbox.zip', files );
