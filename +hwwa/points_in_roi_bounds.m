function [ib, labels, new_to_orig] = points_in_roi_bounds(points, image_rois, point_labels, rois, roi_labels)

assert_ispair( image_rois, point_labels );
assert_ispair( points, point_labels );
assert_ispair( rois, roi_labels );

validateattributes( image_rois, {'double'}, {'2d', 'ncols', 4}, mfilename, 'image rois' );
validateattributes( points, {'double'}, {'2d', 'ncols', 2}, mfilename, 'points' );

image_files = categorical( point_labels, 'image_file' );
[each_roi_labels, roi_I] = keepeach( roi_labels', 'roi' );

ib = cell( size(roi_I) );
labels = cell( size(ib) );

parfor i = 1:numel(roi_I)
  shared_utils.general.progress( i, numel(roi_I) );
  
  roi_ind = roi_I{i};
  roi_image_files = categorical( roi_labels, 'image_file', roi_ind );
  
  tmp_labels = point_labels';
  tmp_ib = false( size(image_files) );
  
  for j = 1:numel(image_files)
    match_ind = roi_image_files == image_files(j);
    
    if ( nnz(match_ind) == 0 )
      continue;
    end
    
    match_roi = rois(roi_ind(match_ind), :);
    
    if ( rows(match_roi) > 1 )
      match_roi = unique( match_roi, 'rows' );
    end
    
    if ( rows(match_roi) ~= 1 )
      error( 'More than one roi matched "%s".', image_files{j} );
    end
    
    x = points(j, 1);
    y = points(j, 2);
    image_roi = image_rois(j, :);
    targ_roi = hwwa.proportional_to_absolute_roi( image_roi, match_roi );
    
    tmp_ib(j) = x >= targ_roi(1) && x <= targ_roi(3) && ...
      y >= targ_roi(2) && y <= targ_roi(4);
  end
  
  labels{i} = join( tmp_labels, prune(each_roi_labels(i)) );
  ib{i} = tmp_ib;
end

ib = vertcat( ib{:} );
labels = vertcat( fcat, labels{:} );
new_to_orig = repmat( columnize(1:numel(image_files)), numel(roi_I), 1 );

assert_ispair( ib, labels );
assert_ispair( new_to_orig, labels );

end