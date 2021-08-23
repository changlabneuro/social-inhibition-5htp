function defaults = edf_trials(varargin)

defaults = hwwa.get_common_make_defaults( varargin{:} );

defaults.look_back = -500;  % ms
defaults.look_ahead = 500;
defaults.event = '';
defaults.event_subdir = 'el_events';
defaults.output_directory = 'edf_trials';

end