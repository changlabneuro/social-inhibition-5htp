function edf_file = edfs(files, conf)

%   EDFS -- Create edf file.
%
%     file = ... edfs( files ) creates the edf intermediate file from the 
%     unified intermediate file contained in `files`. `files` is a
%     containers.Map or struct with an entry for 'unified'.
%
%     file = ... edfs( ..., conf ) uses `conf` to obtain the root data
%     path, instead of the saved config file.
%
%     This function depends on the Edf2Mat utility provided by SR-research.
%
%     IN:
%       - `files` (containers.Map, struct)
%       - `conf` (struct) |OPTIONAL|
%     FILES:
%       - 'unified'
%     OUT:
%       - `edf_file` (struct)

assert( ~isempty(which('Edf2Mat')), 'No Edf2Mat conversion utility found.' );

if ( nargin < 2 || isempty(conf) )
  conf = hwwa.config.load();
else
  hwwa.util.assertions.assert__is_config( conf );
end

hwwa.validatefiles( files, 'unified' );

data_root = hwwa.dataroot( conf );

unified_file = shared_utils.general.get( files, 'unified' );
unified_filename = hwwa.try_get_unified_filename( unified_file );

unified_id = get_unified_id( unified_filename );
raw_subdirectory = unified_file.raw_subdir;

[edf_fullfile, edf_filename] = get_edf_fileparts( data_root, raw_subdirectory, unified_id );

assert( shared_utils.io.fexists(edf_fullfile), 'edf file "%s" does not exist.' ...
  , edf_filename );

edf_obj = Edf2Mat( edf_fullfile );

samples = edf_obj.Samples;
events = edf_obj.Events;

edf_file = struct();

edf_file.unified_filename = unified_filename;
edf_file.Samples = samples;
edf_file.Events = events;

end

function [edf_fullfile, edf_filename] = get_edf_fileparts(data_p, raw_subdirectory, unified_id)

edf_filename = [ unified_id, '.edf' ];
edf_fullfile = fullfile( data_p, raw_subdirectory, edf_filename );

end

function id = get_unified_id(unified_filename)

mat_ind = strfind( unified_filename, '.mat' );
  
assert( ~isempty(mat_ind) && mat_ind > 1 ...
  , 'No ".mat" extension found in "%s".', unified_filename );

id = unified_filename(1:mat_ind-1);

end