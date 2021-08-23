function obj = get_looped_make_runner(params)

%   GET_LOOPED_MAKE_RUNNER -- Get pre-configured LoopedMakeRunner instance.
%
%     obj = hwwa.get_looped_make_runner(); returns a valid LoopedMakeRunner
%     instance, configured in a manner expected by hwwa.make_* functions. 
%
%     Relevant property values are defined by hwwa.get_common_make_defaults.
%
%     obj = hwwa.get_looped_make_runner( PARAMS ) uses struct `PARAMS` to 
%     configure the object, instead of the default values returned by
%     hwwa.get_common_make_defaults.
%
%     obj = hwwa.get_looped_make_runner( PARAMS, 'name1', value1, ... )
%     assigns value1 to field 'name1' of `PARAMS`, and so on, and then uses
%     the assigned values to configure `obj`. If `PARAMS` is an empty array
%     ([]), values are assigned to the defaults returned by
%     hwwa.get_common_make_defaults.
%
%     See also shared_utils.pipeline.LoopedMakeRunner,
%       hwwa.get_common_make_defaults
%
%     IN:
%       - `params` (struct)
%     OUT:
%       - `obj` (shared_utils.pipeline.LoopedMakeRunner)

if ( nargin < 1 || isempty(params) )
  params = hwwa.get_common_make_defaults();
end

if ( nargin > 1 )
  params = hwwa.parsestruct( params, varargin );
end

obj = shared_utils.pipeline.LoopedMakeRunner;

% Whether to attempt to save output.
obj.save = params.save;

% Whether to attempt to run the function in parallel.
obj.is_parallel = params.is_parallel;

% Whether to allow existing files to be overwritten.
obj.overwrite = params.overwrite;

% Whether to keep the output in memory after saving, and return as part of
% `results`.
obj.keep_output = params.keep_output;

% Function that restricts the list of files in a directory to those
% containing string(s). However, if `files_containing` is empty, then all 
% files are used.
obj.filter_files_func = @(x) hwwa.files_containing( x, params.files_containing );

% Function that obtains the unified_filename from a loaded file.
obj.get_identifier_func = @(x, y) hwwa.try_get_unified_filename( x );

% Controls the verbosity of output to the console.
obj.log_level = params.log_level;

% Optionally sets error-handling behavior.
if ( ~strcmp(params.error_handler, 'default') )
  obj.set_error_handler( params.error_handler );
end

% Optionally avoid processing file identifiers that are already present
% in the object's output directory.
if ( params.skip_existing )
  obj.set_skip_existing_files();
end

end