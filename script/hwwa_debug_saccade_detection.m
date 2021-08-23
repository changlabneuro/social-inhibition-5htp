aligned_outs = hwwa_load_edf_aligned( ...
    'start_event_name', 'go_target_onset' ...
  , 'look_back', 75 ...
  , 'look_ahead', 400 ...
  , 'is_parallel', true ...
  , 'files_containing', hwwa.approach_avoid_files() ...
  , 'error_handler', 'error' ...
);

%%

saccade_outs = hwwa_saccade_info( aligned_outs );

%%

nogo_cue_rois = unique( cat_expanded(1, {aligned_outs.rois.nogo_cue}), 'rows' );
go_target_rois = unique( cat_expanded(1, {aligned_outs.rois.go_target}), 'rows' );
go_target_rois = arrayfun( @(x) go_target_rois(x, :), 1:rows(go_target_rois), 'un', 0 );

assert( rows(nogo_cue_rois) == 1, 'Expected only one roi for nogo_cue; got %d', rows(nogo_cue_rois) );
nogo_cue_rois = {nogo_cue_rois};

%  , 'pre_plot_func', 'centroid_relative_distance' ...

hwwa_plot_saccade_info( saccade_outs ...
  , 'do_save', true ...
  , 'scatter_xlims', [0, 1600] ...
  , 'scatter_ylims', [0, 1200] ...
  , 'scatter_rois', go_target_rois ...
  , 'scatter_points_to_right', false ...
);

%%  evaluate

mask = find( aligned_outputs.labels, {'go_trial', 'correct_true', 'tar'} );
i = 1;

while ( true )
  ind = mask(i);
  starts = start_stops{ind}(:, 1);
  stops = start_stops{ind}(:, 2);
  durs = stops - starts;
  
%   starts = starts(durs > 50);
%   stops = stops(durs > 50);
  
  hold off;
  plot( smooth_func(x_deg(ind, :)) );
  hold on;
  plot( smooth_func(y_deg(ind, :)) );
  plot( x_deg(ind, :) );
  plot( y_deg(ind, :) );
  shared_utils.plot.add_vertical_lines( gca, starts, 'k' );
  shared_utils.plot.add_vertical_lines( gca, stops, 'k--' );
  i = i + 1;
  
  z = input( '' );
end