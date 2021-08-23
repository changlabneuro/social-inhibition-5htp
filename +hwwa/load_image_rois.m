function outs = load_image_rois(varargin)

filename = 'ROI Coordinates_one_sheet_2.xlsx';
file_path = fullfile( hwwa.dataroot(varargin{:}), 'roi_coordinates', filename );

num_sets = 4;
rois = { 'face', 'eye', 'mouth' };
num_rois = numel( rois );

rects = cell( num_sets, 1 );
all_components = cell( size(rects) );
labels = cell( size(rects) );
categories = { 'image_file', 'roi' };

for i = 1:num_sets
  set_name = sprintf( 'Set %d- Update', i );
  [~, ~, raw_xls] = xlsread( file_path, set_name );
  header = raw_xls(1, :);
  
  is_x = contains_str( header, ' x' );
  is_y = contains_str( header, ' y' );
  is_w = contains_str( header, ' width' );
  is_h = contains_str( header, ' height' );
  
  assert( nnz(is_x) == num_rois && nnz(is_y) == num_rois && nnz(is_w) && nnz(is_h) ...
    , 'Mismatched components.' );
  
  num_rects = size( raw_xls, 1 ) - 1;
  
  tmp_rects = cell( num_rects * num_rois, 1 );
  tmp_components = cell( size(tmp_rects) );
  tmp_labels = cell( num_rects * num_rois, 2 );
  stp = 1;
  to_keep = true( size(tmp_rects) );
  
  for j = 2:size(raw_xls, 1)
    image_filename = raw_xls{j, 1};
    
    for k = 1:numel(rois)
      if ( ~ischar(image_filename) )
        to_keep(stp) = false;
        stp = stp + 1;
        continue;
      end
      
      if ( ~endsWith(image_filename, '.jpg') )
        image_filename = sprintf( '%s.jpg', image_filename );
      end
      
      roi = rois{k};
      
      is_roi = contains_str( header, roi );
      
      is_x_full = is_x & is_roi;
      is_y_full = is_y & is_roi;
      is_w_full = is_w & is_roi;
      is_h_full = is_h & is_roi;
      
      assert( nnz(is_x_full) == 1 && nnz(is_y_full) == 1 && nnz(is_w_full) == 1 && nnz(is_h_full) == 1 ...
        , 'Mismatched components.' );
      
      x = raw_xls{j, is_x_full};
      y = raw_xls{j, is_y_full};
      w = raw_xls{j, is_w_full};
      h = raw_xls{j, is_h_full};
      
      rect = [ x, y, x+w, y+h ];
      components = [ x, y, w, h ];

      assert( ~isempty(rect) && ~isempty(components) );
      
      tmp_rects{stp} = rect;
      tmp_components{stp} = components;
      tmp_labels(stp, :) = { image_filename, roi };
      
      stp = stp + 1;
    end
  end
  
  all_components{i} = vertcat( tmp_components{to_keep} );
  rects{i} = vertcat( tmp_rects{to_keep} );
  labels{i} = tmp_labels(to_keep, :);
end

rects = vertcat( rects{:} );
all_components = vertcat( all_components{:} );
labels = vertcat( labels{:} );

assert( rows(rects) == rows(all_components) );
assert( rows(rects) == rows(labels) );

outs = struct();
outs.rects = rects;
outs.components = all_components;
outs.labels = labels;
outs.categories = categories;

end

function tf = contains_str(header, s)

tf = cellfun( @(x) ischar(x) && ~isempty(strfind(lower(x), s)), header );

end