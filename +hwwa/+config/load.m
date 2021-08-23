
function out = load(cmd)

%   LOAD -- Load the config file.
%
%     out = hwwa.config.load() loads the current config file.
%     out = hwwa.config.load( '-clear' ) clears the cache of the config
%     file before loading. This isn't ever necessary unless manually
%     updating the `config.mat` file via the filesystem.
%
%     OUT:
%       - `out` (struct)

% cache `conf`, unless it is updated via a call to `hwwa.config.save`
persistent conf;

if ( nargin > 0 )
  assert( strcmp(cmd, '-clear'), 'Invalid command. Valid options are: -clear.' );  
  conf = [];
end

if ( isempty(conf) )
  conf = do_load();
end

out = conf;

end

function conf = do_load()

const = hwwa.config.constants();

config_folder = const.config_folder;
config_file = const.config_filename;

config_filepath = fullfile( config_folder, config_file );

if ( exist(config_filepath, 'file') ~= 2 )
  fprintf( '\n Creating config file ...' );
  conf = hwwa.config.create();
  hwwa.config.save( conf );
end

conf = load( config_filepath );
conf = conf.(char(fieldnames(conf)));

try
  hwwa.util.assertions.assert__is_config( conf );
catch err
  fprintf( ['\n Could not load the config file saved at ''%s'' because ' ...
    , ' it is not recognized as a valid config file.\n\n'], config_filepath );
  throwAsCaller( err );
end

end