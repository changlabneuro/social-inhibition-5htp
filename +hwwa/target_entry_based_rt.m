function rts = target_entry_based_rt(aligned_outs, aligned_out_mask)

if ( nargin < 2 )
  aligned_out_mask = rowmask( aligned_outs.labels );
end

padded_go_targ_rois = hwwa.linearize_roi( aligned_outs, 'go_target_padded' );
go_targ_displacement = hwwa.linearize_roi( aligned_outs, 'go_target_displacement' );

rois = padded_go_targ_rois;
displacement = go_targ_displacement;

directional_go_targ_rois = ...
  hwwa.direction_dependent_go_target_rois( rois, aligned_outs.labels ...
  , shared_utils.rect.width(displacement), shared_utils.rect.height(displacement) );

t = aligned_outs.params.look_back:aligned_outs.params.look_ahead;
assert( numel(t) == size(aligned_outs.x, 2), 'Time does not match position' );

use_time_of_first_entry_in_target = true;

rts = nan( rows(aligned_outs.labels), 1 );

for i = 1:numel(aligned_out_mask)
  ind = aligned_out_mask(i);
  roi = directional_go_targ_rois(ind, :);
  
  aligned_x = aligned_outs.x(ind, :);
  aligned_y = aligned_outs.y(ind, :);
  all_ib = arrayfun( @(x, y) bfw.bounds.rect(x, y, roi), aligned_x, aligned_y );
 
  [starts, durs] = shared_utils.logical.find_islands( all_ib );

  if ( ~isempty(starts) )
    if ( use_time_of_first_entry_in_target )
      start_ind = starts(1);
    else
      start_ind = starts(1) + durs(1) - 1;
    end

    rts(ind) = t(start_ind);
  end
end

end