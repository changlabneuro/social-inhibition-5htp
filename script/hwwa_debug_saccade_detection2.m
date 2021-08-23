conf = hwwa.config.load();
aligned_outs = hwwa_approach_avoid_go_target_on_aligned( 'config', conf );
roi_info = hwwa.load_image_rois();

%%

saccade_outs = hwwa_saccade_info( aligned_outs );

%%

go_target_rois = hwwa.linearize_roi( aligned_outs, 'go_target' );
nogo_cue_rois = hwwa.linearize_roi( aligned_outs, 'nogo_cue' );
padded_go_targ_rois = hwwa.linearize_roi( aligned_outs, 'go_target_padded' );
go_targ_displacement = hwwa.linearize_roi( aligned_outs, 'go_target_displacement' );

saccade_rois = go_target_rois(saccade_outs.aligned_to_saccade_ind, :);

%%

rects = hwwa.proportional_rois( roi_info.rects );
roi_labels = fcat.from( roi_info );
rects = hwwa.add_full_image_proportional_roi( rects, roi_labels );

[new_rects, new_labs] = hwwa.add_top_bottom_half_face_rois( rects, roi_labels );
roi_labels = [ roi_labels'; new_labs ];
rects = [ rects; new_rects ];

[ib, labels] = hwwa.points_in_roi_bounds( saccade_outs.saccade_stop_points ...
  , saccade_rois, saccade_outs.labels, rects, roi_labels );

saccade_outs.in_bounds_roi_labels = labels';
saccade_outs.in_bounds_roi = ib;

%%

do_save = true;

norm_combs = [ true ];
per_monk_combs = [ false ];
per_image_cat_combs = [ false ];
per_correct_types = false;
comb_inds = dsp3.numel_combvec( norm_combs, per_monk_combs ...
  , per_image_cat_combs, per_correct_types );

for i = 1:size(comb_inds, 2)
  inputs = {};
  
  is_norm = norm_combs(comb_inds(1, i));
  is_per_monk = per_monk_combs(comb_inds(2, i));
  is_per_image_cat = per_image_cat_combs(comb_inds(3, i));
  is_per_correct = per_correct_types(comb_inds(4, i));
  
  base_subdir = '';
  base_subdir = sprintf( '%s-%s', base_subdir, ternary(is_norm, 'norm', 'non-norm') );
  base_subdir = sprintf( '%s%s%s', base_subdir, filesep, ternary(is_per_monk, 'per-monk', 'across-monks') );
  
  if ( is_norm )
    inputs = [ inputs, {'pre_plot_func', 'saline_normalize'} ];
  end
  
%   funcs = {'in_bounds_roi_proportions', 'saccade_landing_point_stddev'};
  funcs = {'in_bounds_roi_proportions'};
%   funcs = { 'saccade_landing_point_scatter' };
  mask_func = @(labels, mask) fcat.mask(labels, mask ...
    , @find, {'face_bot_half', 'face_top_half'} ...
  );

  hwwa_plot_saccade_info( saccade_outs ...
    , 'do_save', true ...
    , 'scatter_xlims', [0, 1600] ...
    , 'scatter_ylims', [0, 1200] ...
    , 'scatter_rois', go_target_rois ...
    , 'scatter_points_to_right', false ...
    , 'normalize', is_norm ...
    , 'funcs', funcs ...
    , 'base_subdir', base_subdir ...
    , 'per_monkey', is_per_monk ...
    , 'per_image_category', is_per_image_cat ...
    , 'prefer_boxplot', true ...
    , 'per_correct_type', is_per_correct ...
    , 'mask_func', mask_func ...
    , inputs{:} ...
  );
end


%%

is_proportional_rects = true;

clf();

use_rects = ternary( is_proportional_rects, rects, roi_info.rects );

image_filename = cellstr( roi_labels, 'image_file', 1 );
image_mask = find( roi_labels, image_filename );
% image_mask = find( roi_labels, {'eye', 'face'}, image_mask );

roi_ind = findall( roi_labels, 'roi', image_mask );

image_roi = saccade_rois(1, :);
colors = spring( numel(roi_ind) );

for i = 1:numel(roi_ind)
  rect = use_rects(roi_ind{i}, :);
  
  if ( is_proportional_rects )
    rect = hwwa.proportional_to_absolute_roi( image_roi, rect );
  end
  
  hs = bfw.plot_rect_as_lines( gca, rect );
  set( hs, 'color', colors(i, :) );
end

%%  

