function conf = create(do_save)

%   CREATE -- Create the config file. 
%
%     Define editable properties of the config file here.
%
%     IN:
%       - `do_save` (logical) -- Indicate whether to save the created
%         config file. Default is `false`

if ( nargin < 1 ), do_save = false; end

const = hwwa.config.constants();

conf = struct();

%   ID
conf.(const.config_id) = true;

project_dir = hwwa.util.get_project_folder();

%   PATHS
PATHS = struct();
PATHS.data_root = fullfile( project_dir, 'data' );
PATHS.repositories = fileparts( project_dir );
PATHS.raw_subdirectory = 'raw';

%   DEPENDENCIES
DEPENDS = struct();
DEPENDS.repositories = { 'shared_utils', 'eyelink', 'plexon' ...
  , 'spike_helpers', 'ms_run', 'dsp3', 'categorical/api/matlab' };
DEPENDS.others = { '' };

%   PLEX
PLEX = struct();
PLEX.sync_channel = 'AI01';

%   EXPORT
conf.PATHS = PATHS;
conf.DEPENDS = DEPENDS;
conf.PLEX = PLEX;

if ( do_save )
  hwwa.config.save( conf );
end

end