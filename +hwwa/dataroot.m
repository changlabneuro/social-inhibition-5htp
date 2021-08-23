function p = dataroot(conf)

%   DATAROOT -- Get the absolute path to the root data directory.
%
%     p = hwwa.dataroot() returns the path to the root data directory, as
%     defined in the saved config file.
%
%     p = hwwa.dataroot( CONF ) uses the config file `CONF`, instead of the
%     saved config file.
%
%     IN:
%       - `conf` (struct) |OPTIONAL|
%     OUT:
%       - `p` (char)

if ( nargin < 1 || isempty(conf) )
  conf = hwwa.config.load();
else
  hwwa.util.assertions.assert__is_config( conf );
end

p = conf.PATHS.data_root;

end