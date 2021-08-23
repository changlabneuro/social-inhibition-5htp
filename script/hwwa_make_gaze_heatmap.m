function [heat_maps, heat_map_labs, x_edges, y_edges] = ...
  hwwa_make_gaze_heatmap(x, y, labels, heat_map_each, x_lims, y_lims, stp, win, varargin)

defaults = struct();
defaults.mask = rowmask( labels );
params = hwwa.parsestruct( defaults, varargin );

assert_ispair( x, labels );
assert_ispair( y, labels );

mask = params.mask;

heat_map_I = findall( labels, heat_map_each, mask );

heat_maps = {};
heat_map_labs = {};
x_edges = {};
y_edges = {};

parfor i = 1:numel(heat_map_I)  
  shared_utils.general.progress( i, numel(heat_map_I) );
  
  [tmp_map, x_edge, y_edge] = hwwa_gaze_heatmap( x, y, x_lims, y_lims, stp, win, heat_map_I{i} );
  tmp_map = reshape( tmp_map, [1, size(tmp_map)] );
  
  heat_map_labs{i} = append1( fcat(), labels, heat_map_I{i} );
  heat_maps{i} = tmp_map;
  x_edges{i} = x_edge;
  y_edges{i} = y_edge;
end

heat_map_labs = vertcat( fcat(), heat_map_labs{:} );
heat_maps = vertcat( heat_maps{:} );

if ( ~isempty(x_edges) )
  x_edges = x_edges{1};
  y_edges = y_edges{1};
else
  x_edges = [];
  y_edges = [];
end

end