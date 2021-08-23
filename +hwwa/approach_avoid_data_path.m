function p = approach_avoid_data_path(plot_params, data_type, varargin)

p = fullfile( hwwa.dataroot(plot_params.config), data_type, 'approach_avoid' ...
  , 'behavior', hwwa.datedir, plot_params.base_subdir, varargin{:} );

end