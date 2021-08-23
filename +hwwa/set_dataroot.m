function [old_p, conf] = set_dataroot(p, conf)

%   SET_DATAROOT -- Set root data directory.
%
%     hwwa.set_dataroot( PATH ) sets the root data directory to the
%     absolute path given by the character vector `PATH`. The current 
%     config file is modified and saved with the new root directory.
%
%     hwwa.set_dataroot( ..., conf ) uses the config file `conf` instead of
%     the current saved config file.
%
%     old_p = hwwa.set_dataroot( ... ) returns the previous data root
%     directory.
%
%     [.., conf] = hwwa.set_dataroot( ... ) also returns the *updated*
%     config file.
%
%     IN:
%       - `path` (char)
%       - `conf` (struct) |OPTIONAL|
%     OUT:
%       - `old_path` (char)
%       - `conf` (struct)

if ( nargin < 2 || isempty(conf) )
  conf = hwwa.config.load();
else
  hwwa.util.assertions.assert__is_config( conf );
end

validateattributes( p, {'char'}, {'scalartext'}, mfilename, 'p' );

old_p = conf.PATHS.data_root;
conf.PATHS.data_root = p;

hwwa.config.save( conf );

end