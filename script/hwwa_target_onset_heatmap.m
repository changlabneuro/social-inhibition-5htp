aligned_outputs = hwwa_load_edf_aligned( ...
    'start_event_name', 'go_target_onset' ...
  , 'look_back', 300 ...
  , 'look_ahead', 750 ...
  , 'is_parallel', true ...
  , 'files_containing', hwwa.approach_avoid_files() ...
);

%%

m = hwwa.get_approach_avoid_mask( aligned_outputs.labels );
[count_labs, each_I] = keepeach( aligned_outputs.labels', {'unified_filename'}, m );
counts = zeros( numel(each_I), 2 );

for i = 1:numel(each_I)
  num_l = numel( find(aligned_outputs.labels, 'center-left', each_I{i}) );
  num_r = numel( find(aligned_outputs.labels, 'center-right', each_I{i}) );
  counts(i, :) = [ num_l, num_r ];
end

count_labs = repset( count_labs, 'target_placement', {'center-left', 'center-right'} );
counts1 = columnize( counts );

pl = plotlabeled.make_common();
axs = pl.bar( counts1, count_labs, 'target_placement', {}, {} );

%%

aligned_outputs = shared_utils.io.fload( ...
  fullfile(hwwa.gid('processed'), 'behavior/heatmap/go_targ_aligned.mat') );

%%

norm_x = aligned_outputs.x / monitor_width_px;
norm_y = 1 - aligned_outputs.y / monitor_height_px;
labels = aligned_outputs.labels';

x_lims = [ -0.2, 1.2 ];
y_lims = [ -0.2, 1.2 ];

stp = 0.01;
win = 0.01;

% mask = fcat.mask( labels, hwwa.get_approach_avoid_mask(labels) ...
%   , @find, {'nogo_trial', 'correct_true'} ...
% );

% mask = fcat.mask( labels, hwwa.get_approach_avoid_mask(labels) ...
%   , @find, {'correct_true'} ...
% );

mask = fcat.mask( labels, hwwa.get_approach_avoid_mask(labels) );

% heat_map_each = { 'trial_type', 'target_image_category', 'scrambled_type', 'monkey', 'drug', 'day' };
% heat_map_each = { 'trial_type', 'monkey', 'drug' };
heat_map_each = { 'trial_type', 'drug', 'correct', 'unified_filename' };

[heat_maps, heat_map_labs, x_edges, y_edges] = hwwa_make_gaze_heatmap( norm_x, norm_y, labels, heat_map_each, x_lims, y_lims, stp, win ...
  , 'mask', mask ...
);

%%

norm_roi = make_normalized_roi( aligned_outputs, 'nogo_cue' );
targ_rois = get_unique_rois( aligned_outputs, 'go_target' );
targ_rois = targ_rois([2, 4], :);

targ_rois = arrayfun( @(x) normalize_roi(targ_rois(x, :)), 1:size(targ_rois, 1), 'un', 0 );

%%  compare 5-htp and saline
  
% test_rois = { {targ_rois{1}, targ_rois{2}} };
% test_roi_names = { 'left-right-target' };

use_heat_map = heat_maps;
use_x_edges = x_edges;
use_y_edges = y_edges;

do_collapse = false;

if ( do_collapse )
  [use_heat_map, use_x_edges, use_y_edges] = collapse_to_right( use_heat_map, use_x_edges, use_y_edges );
end

test_rois = { targ_rois(1), targ_rois(2), {norm_roi} };
test_roi_names = { 'left-target', 'right-target', 'cue' };

rs_tables = [];
rs_labels = fcat();

store_dats = {};
store_labs = {};

curr_trial_type = 'go_trial';

for i = 1:numel(test_rois)
  test_labs = fcat();
  test_dat = [];
  
  for j = 1:numel(test_rois{i})
    collapsed = roi_collapse( use_heat_map, test_rois{i}{j}, use_x_edges, use_y_edges );
    test_mask = find( heat_map_labs, 'correct_true' );
    test_mask = find( heat_map_labs, curr_trial_type, test_mask );
    
    sub_collapsed = collapsed(test_mask); 
    
    append( test_labs, heat_map_labs, test_mask );
    test_dat = [ test_dat; sub_collapsed ];
  end
  
  store_dats{end+1, 1} = test_dat;
  store_labs{end+1, 1} = test_labs';
    
  rs_outs = dsp3.ranksum( test_dat, test_labs', {'trial_type'}, '5-htp', 'saline' );
  
  addsetcat( rs_outs.rs_labels, 'target_name', test_roi_names{i} );
  append( rs_labels, rs_outs.rs_labels );
  rs_tables = [ rs_tables; rs_outs.rs_tables ];
end

assert_ispair( rs_tables, rs_labels );

rs_tables = vertcat( rs_tables{:} );
sig_ind = rs_tables.p < 0.05;
rs_labels(find(sig_ind), {'trial_type', 'target_name'})

if ( true )
  for i = 1:numel(store_dats)
    data_tbl = array2table( store_dats{i}, 'VariableNames', {'fix_pattern_freq'} );
    label_tbl = array2table( cellstr(store_labs{i}), 'VariableNames', getcats(store_labs{i}) );
    save_table_pair( data_tbl, label_tbl, sprintf('fix_pattern_freq_%d_%s', i, curr_trial_type), '' );
  end
end

%%

select_x = x_edges >= 0.2 & x_edges <= 0.75;
% select_y = y_edges >= 0.3 & y_edges <= 0.7;
select_y = y_edges <= 0.72;

select_x = true( size(x_edges) );
select_y = true( size(y_edges) );

subset_dat = heat_maps(:, select_y, select_x);

soc_minus_scr = ...
  @(data, labels, spec) hwwa.social_minus_scrambled( data, labels', setdiff(spec, 'scrambled_type') );

fhtp_minus_sal = ...
  @(data, labels, spec) hwwa.fhtp_minus_saline( data, labels', setdiff(spec, 'drug') );

noop = @(data, labels, spec) deal(data, labels);

max_norm = @hwwa.max_normalize;

ys = y_edges(select_y);
xs = x_edges(select_x);

ys_deg = normalized_edges_to_degrees( ys - 0.5, monitor_height_px );
xs_deg = normalized_edges_to_degrees( xs - 0.5, monitor_width_px );

norm_rois = [ norm_roi, targ_rois ];

for i = 1:numel(norm_rois)
  max_vs = [ monitor_width_px, monitor_height_px, monitor_width_px, monitor_height_px ];
  
  for j = 1:numel(max_vs)
    norm_rois{i}(j) = normalized_edges_to_degrees( norm_rois{i}(j) - 0.5, max_vs(j) );
  end
end

hwwa_plot_target_onset_heatmap( subset_dat, heat_map_labs', ys_deg, xs_deg ...
  , 'do_save', true ...
  , 'before_plot_func', max_norm ...
  , 'base_subdir', 'max-norm' ...
  , 'match_c_lims', true ...
  , 'overlay_rects', norm_rois ...
  , 'c_lims', [0.0, 0.5] ...
);

function w = monitor_width_px()
w = 1600;
end

function h = monitor_height_px()
h = 900;
end

function degs = normalized_edges_to_degrees(edges, max_v)

pxs = edges * max_v;
degs_old = hwwa.run_px2deg( pxs, [] );
degs = hwwa.convert_from_old_degrees_to_corrected_degrees( degs_old );

end

function rois = get_unique_rois(aligned_outputs, roi_name)

rois = hwwa.linearize_roi( aligned_outputs, roi_name );
rois = unique( rois, 'rows' );

end

function norm_roi = normalize_roi(roi)

norm_roi = roi;
norm_roi([1, 3]) = norm_roi([1, 3]) / monitor_width_px;
norm_roi([2, 4]) = norm_roi([2, 4]) / monitor_height_px;

end

function norm_roi = make_normalized_roi(aligned_outputs, roi_name)

unique_rois = get_unique_rois( aligned_outputs, roi_name );
assert( rows(unique_rois) == 1 );
norm_roi = normalize_roi( unique_rois );

end

function m = roi_collapse(heat_map, roi, x, y)

assert( numel(x) == size(heat_map, 3) && numel(y) == size(heat_map, 2) );
x_ind = x >= roi(1) & x < roi(3);
y_ind = y >= roi(2) & y < roi(4);
window = heat_map(:, y_ind, x_ind);
m = sum( sum(window, 3), 2 );
n = sum( x_ind ) * sum( y_ind );

m = m / n;

end

function [h, x, y] = collapse_to_right(heat_map, x, y)

assert( numel(x) == size(heat_map, 3) && numel(y) == size(heat_map, 2) );
ind_r = x >= 0.5 & x < 1;
ind_l = x >= 0 & x < 0.5;

assert( sum(ind_r) == sum(ind_l) );

h_r = heat_map(:, :, ind_r);
h_l = flip( heat_map(:, :, ind_l), 3 );

h = h_r + h_l;
x = x(ind_r);

end

function save_table_pair(data_tbl, label_tbl, subdir, fname)

assert( size(data_tbl, 1) == size(label_tbl, 1) );

params = hwwa.get_common_make_defaults( hwwa.get_common_plot_defaults() );
save_p = fullfile( hwwa.approach_avoid_data_path(params, 'export'), subdir );

shared_utils.io.require_dir( save_p );
dsp3.writetable( data_tbl, fullfile(save_p, sprintf('data__%s.csv', fname)) );
dsp3.writetable( label_tbl, fullfile(save_p, sprintf('labels__%s.csv', fname)) );

end
