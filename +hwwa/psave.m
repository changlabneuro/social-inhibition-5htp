function psave( pathstr, var, varargin )

%   PSAVE -- Save variable in parfor loop.
%
%     psave( 'eg.mat', 10 ) writes 10 to the file 'eg.mat'. In the .mat
%     file, the variable is bound to the name 'var'.
%
%     psave( 'eg.mat', 10, 'X' ) works as above, but binds the variable to
%     the name 'X'. In this way, load('eg.mat') will import a variable
%     named 'X' to the workspace.
%
%     IN:
%       - `pathstr` (char)
%       - `var` (/any/)
%       - `alias` (char) |OPTIONAL|

if ( nargin < 3 )
  save( pathstr, 'var', varargin{:} );
else
  alias = varargin{1};
  varargin(1) = [];
  eval( sprintf('%s=var;', alias) );
  save( pathstr, alias, varargin{:} );
end

end