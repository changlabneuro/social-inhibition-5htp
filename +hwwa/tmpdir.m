function p = tmpdir(varargin)

%   TMPDIR -- Temporary directory path.
%
%     See also hwwa.dataroot

dr = hwwa.dataroot( varargin{:} );
p = fullfile( dr, 'tmp' );

end