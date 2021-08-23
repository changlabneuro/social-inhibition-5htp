function rois = linearize_roi(aligned_outs, roi_name)

base_rois = cat_expanded( 1, columnize({aligned_outs.rois.(roi_name)}) );
rois = base_rois(aligned_outs.roi_indices, :);

end