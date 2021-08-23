function p = get_intermediate_dir(kind, conf)

%   GET_INTERMEDIATE_DIR -- Get the absolute path to an intermediate directory.
%
%     p = ... get_intermediate_dir( KIND ); returns the absolute path to
%     the intermediate directory `KIND`, using the root data directory
%     given by the saved config file.
%
%     p = ... get_intermediate_dir( ..., CONF ); uses the config file
%     `CONF` to generate the absolute path, instead of the saved config
%     file.
%
%     See also hwwa.load_one_intermediate
%
%     IN:
%       - `kind` (cell array of strings, char)
%       - `conf` (struct) |OPTIONAL|
%     OUT:
%       - `P` (char)

if ( nargin < 1 ), kind = ''; end
if ( nargin < 2 || isempty(conf) )
  conf = hwwa.config.load();
else
  hwwa.util.assertions.assert__is_config( conf );
end

p = shared_utils.io.fullfiles( conf.PATHS.data_root, 'intermediates', kind );

if ( ischar(kind) )
  p = char( p );
end

end