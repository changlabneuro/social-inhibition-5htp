function unified_file = unified(files, unified_filename, conf)

%   UNIFIED -- Create unified file.
%
%     file = ... unified(files, unified_filename) creates the unified 
%     intermediate file from the raw data file contained in `files`, 
%     assigning that file the char identifier `unified_filename`. `files`
%     is a struct or containers.Map with an entry for 'raw', the name
%     of the subdirectory in which raw data files were presumed to be saved.
%
%     file = ... unified( ..., conf ) uses the config file `conf` to get
%     the name of the subdirectory in which raw-data is stored, instead of
%     the saved config file. The value of this `raw_subdirectory` field
%     must be consistent with the entry into the `files` aggregate.
%
%     IN:
%       - `files` (containers.Map, struct)
%       - `unified_filename` (char)
%       - `conf` (struct) |OPTIONAL|
%     FILES:
%       - <raw>
%     OUT:
%       - `unified_file` (struct)

if ( nargin < 3 || isempty(conf) )
  conf = hwwa.config.load();
else
  hwwa.util.assertions.assert__is_config( conf );
end

raw_subdir = conf.PATHS.raw_subdirectory;

% Ensure file associated with `raw_subdir` is present.
hwwa.validatefiles( files, raw_subdir );

input_p = fullfile( hwwa.dataroot(conf), raw_subdir );

unified_file = shared_utils.general.get( files, raw_subdir );
unified_id = shared_utils.io.filenames( unified_filename );

raw_data = errors_to_string( unified_file.DATA );
  
reg_map = get_region_map( input_p, unified_id );

unified_file.DATA = raw_data;
unified_file.unified_filename = unified_filename;
unified_file.unified_id = unified_id;
unified_file.raw_subdir = raw_subdir;
unified_file.region_map = reg_map;

end

function s = errors_to_string(s)

fields = fieldnames( s );

error_prefix = 'error__';

errs = find( cellfun(@(x) ~isempty(strfind(x, error_prefix)), fields) );

for idx = 1:numel(s)

  str = '';

  for i = 1:numel(errs)
    f = fields{errs(i)};
    val = s(idx).(f);

    if ( val )
      assert( isempty(str), 'More than one error for this trial.' );
      str = f(numel(error_prefix)+1:end);
    end
  end

  if ( isempty(str) )
    str = 'no_errors';
  end

  s(idx).error = str;
end

s = rmfield( s, fields(errs) );

end

function reg_map = get_region_map(input_p, unified_id)

reg_map = [];

reg_file = fullfile( input_p, [unified_id, '.regions'] );

if ( ~shared_utils.io.fexists(reg_file) )
  return;
end

regs = jsondecode( fileread(reg_file) );

reg_map = structfun( @shared_utils.general.json_channels2num, regs, 'un', false );

end