function defaults = get_common_make_defaults(assign_to)

%   GET_COMMON_MAKE_DEFAULTS -- Get common default values for hwwa.make_* 
%     functions.
%
%     IN:
%       - `assign_to` (struct) |OPTIONAL|
%     OUT:
%       - `defaults` (struct)

if ( nargin == 0 )
  defaults = struct();
else
  defaults = assign_to;
end

defaults.loop_runner = [];
defaults.files = [];
defaults.files_containing = [];
defaults.overwrite = false;
defaults.append = true;
defaults.save = true;
defaults.log_level = 'info';
defaults.is_parallel = true;
defaults.keep_output = false;
defaults.error_handler = 'default';
defaults.skip_existing = false;
defaults.config = hwwa.config.load();

end