function [rois, roi_labels] = add_full_image_proportional_roi(rois, roi_labels)

assert_ispair( rois, roi_labels );

roi_I = findall( roi_labels, 'roi' );

if ( ~isempty(roi_I) )
  first_roi_ind = roi_I{1};
  orig_num_rows = rows( roi_labels );
  
  append( roi_labels, roi_labels, first_roi_ind );
  setcat( roi_labels, 'roi', 'image', (orig_num_rows+1):rows(roi_labels) );
  
  new_rois = repmat( [0, 0, 1, 1], numel(first_roi_ind), 1 );
  rois = [ rois; new_rois ];
end

end