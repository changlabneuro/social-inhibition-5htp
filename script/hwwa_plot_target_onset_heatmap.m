function hwwa_plot_target_onset_heatmap(heat_maps, heat_map_labs, x, y, varargin)

defaults = hwwa.get_common_plot_defaults( hwwa.get_common_make_defaults() );
defaults.before_plot_func = @(varargin) deal(varargin{1:nargout});
defaults.match_c_lims = true;
defaults.overlay_rects = [];
defaults.mask_func = @hwwa.default_mask_func;
defaults.c_lims = [];
defaults.custom = false;
defaults.pcats = {};
defaults.fcats = {};
defaults.subdir = '';

params = hwwa.parsestruct( defaults, varargin );
mask = params.mask_func( heat_map_labs, rowmask(heat_map_labs) );

if ( params.custom )
  custom( heat_maps(mask, :, :), prune(heat_map_labs(mask)), x, y, params );
else
  % per_image_cat( heat_maps, heat_map_labs', x, y, params );
  across_image_cat( heat_maps(mask, :, :), prune(heat_map_labs(mask)), x, y, params );
end

end

function custom(heat_maps, heat_map_labs, x, y, params)

plot_heat_map( heat_maps, heat_map_labs, x, y ...
  , params.fcats, params.pcats, params.subdir, params );

end

function across_image_cat(heat_maps, heat_map_labs, x, y, params)

assert_ispair( heat_maps, heat_map_labs );

for i = 1
  subdir = 'across_image_category';
  
  if ( i == 2 )
    collapsecat( heat_map_labs, 'monkey' );
    subdir = sprintf( '%s-across_monkeys', subdir );
  else
    subdir = sprintf( '%s-per_monkey', subdir );
  end
  
  fcats = { 'monkey', 'correct' };
  pcats = { 'trial_type', 'drug', 'monkey', 'correct' };

  plot_heat_map( heat_maps, heat_map_labs, x, y, fcats, pcats, subdir, params );
end

end

function per_image_cat(heat_maps, heat_map_labs, x, y, params)

for i = 1:2
  subdir = 'per_image_category';
  
  if ( i == 2 )
    collapsecat( heat_map_labs, 'monkey' );
    subdir = sprintf( '%s-across_monkeys', subdir );
  else
    subdir = sprintf( '%s-per_monkey', subdir );
  end
  
  fcats = { 'monkey', 'target_image_category' };
  pcats = { 'trial_type', 'drug', 'monkey', 'target_image_category', 'scrambled_type' };

  plot_heat_map( heat_maps, heat_map_labs, x, y, fcats, pcats, subdir, params );
end

end

function r = rect_to_spectrogram_rect(rect, x, y)

nearest_x0 = nearest_bin( rect(1), x );
nearest_x1 = nearest_bin( rect(3), x );
nearest_y0 = nearest_bin( rect(2), y );
nearest_y1 = nearest_bin( rect(4), y );

x0 = to_bin_with_fraction( x, nearest_x0, rect(1) );
x1 = to_bin_with_fraction( x, nearest_x1, rect(3) );
y0 = to_bin_with_fraction( y, nearest_y0, rect(2) );
y1 = to_bin_with_fraction( y, nearest_y1, rect(4) );

r = [ x0, y0, x1, y1 ];

end

function r = transpose_rect(rect)

r = [ rect(2), rect(1), rect(4), rect(3) ];

end

function r = flip_rect_ud(rect, max_y)

r = [ rect(1), max_y - rect(2), rect(3), max_y - rect(4) ];

end

function adjusted_component = to_bin_with_fraction(bins, nearest_ind, component)

frac = (component - bins(nearest_ind)) / (bins(nearest_ind+1) - bins(nearest_ind));
adjusted_component = frac + nearest_ind;

end

function ind = nearest_bin(component, bins)

ind = find( bins > component, 1 ) - 1;

end

function plot_heat_map(heat_maps, heat_map_labs, x, y, fcats, pcats, subdir, params)

conf = params.config;
do_save = params.do_save;

spec = unique( cshorzcat(fcats, pcats) );

fig_I = findall_or_one( heat_map_labs, fcats );

for i = 1:numel(fig_I)
  pl = plotlabeled.make_spectrogram( x, y );
  pl.smooth_func = @(x) imgaussfilt(x, 2);
  pl.add_smoothing = true;
  pl.match_c_lims = params.match_c_lims;
  
  if ( numel(fig_I{i}) == 2 )
    pl.shape = [1, 2];
  end

  plt = heat_maps(fig_I{i}, :, :);
  plt_labs = prune( heat_map_labs(fig_I{i}) );
  
  [plt, plt_labs] = params.before_plot_func( plt, plt_labs, spec );
  
  axs = pl.imagesc( plt, plt_labs, pcats );
  
  if ( ~isempty(params.overlay_rects) )
    for j = 1:numel(params.overlay_rects)
      rect = rect_to_spectrogram_rect( params.overlay_rects{j}, y, x );
      rect = flip_rect_ud( rect, numel(x) );
      hs = arrayfun( @(x) bfw.plot_rect_as_lines(x, rect), axs, 'un', 0 );
      cellfun( @(x) set(x, 'color', zeros(1, 3)), hs );
      cellfun( @(x) set(x, 'linewidth', 2), hs );
    end
  end
  
  if ( ~isempty(params.c_lims) )
    shared_utils.plot.set_clims( axs, params.c_lims );
  end
  
  shared_utils.plot.tseries_xticks( axs, round(y), 10 );
  shared_utils.plot.fseries_yticks( axs, round(x), 10 );
  
  if ( do_save )
    save_p = fullfile( hwwa.dataroot(conf), 'plots', 'approach_avoid' ...
      , 'behavior', hwwa.datedir, 'heat_map', params.base_subdir, subdir );
    
    shared_utils.plot.fullscreen( gcf );
    dsp3.req_savefig( gcf, save_p, plt_labs, [fcats, pcats] );
  end
end

end