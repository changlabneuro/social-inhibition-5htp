function [figs, all_axs, fig_I] = plot_grouped_histogram(data, labels, fcats, pcats, gcats, varargin)

assert_ispair( data, labels );
validateattributes( data, {'double'}, {'column'}, mfilename, 'data' );

defaults = struct();
defaults.hist_inputs = {};
defaults.mask_func = @hwwa.default_mask_func;

params = hwwa.parsestruct( defaults, varargin );
mask = hwwa.make_mask( labels, params.mask_func );

fig_I = findall_or_one( labels, fcats, mask );
figs = cell( numel(fig_I), 1 );
all_axs = cell( size(figs) );

if ( isempty(pcats) )
  pcats = fcats;
end

for i = 1:numel(fig_I)
  fig = figure(i);
  figs{i} = fig;
  fig_ind = fig_I{i};
  all_axs{i} = gobjects( 0 );
  
  [panel_I, panel_C] = findall( labels, pcats, fig_ind );
  sub_shape = plotlabeled.get_subplot_shape( numel(panel_I) );
  
  for j = 1:numel(panel_I)
    ax = subplot( sub_shape(1), sub_shape(2), j );
    cla( ax );
    hold( ax, 'on' );
    
    [group_I, group_C] = findall( labels, gcats, panel_I{j} );
    hs = gobjects( numel(group_I), 1 );    
    for k = 1:numel(group_I)
      hs(k) = histogram( ax, data(group_I{k}), params.hist_inputs{:} );
    end
    
    legend( hs, strrep(fcat.strjoin(group_C), '_', ' ') );
    title( ax, strrep(fcat.strjoin(panel_C), '_', ' ') );
    all_axs{i}(end+1, 1) = ax;
  end
end

end