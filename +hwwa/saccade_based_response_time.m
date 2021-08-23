function rts = saccade_based_response_time(aligned_outs, saccade_outs, saccade_out_mask, varargin)

defaults = struct();
defaults.use_time_of_first_entry_in_target = false;
defaults.allow_out_of_bounds_end_point = true;

params = hwwa.parsestruct( defaults, varargin );

if ( nargin < 3 )
  saccade_out_mask = rowmask( saccade_outs.labels );
end

padded_go_targ_rois = hwwa.linearize_roi( aligned_outs, 'go_target_padded' );
go_targ_displacement = hwwa.linearize_roi( aligned_outs, 'go_target_displacement' );

rois = padded_go_targ_rois(saccade_outs.aligned_to_saccade_ind, :);
displacement = go_targ_displacement(saccade_outs.aligned_to_saccade_ind, :);
stops = saccade_outs.saccade_stop_points;

directional_go_targ_rois = ...
  hwwa.direction_dependent_go_target_rois( rois, saccade_outs.labels ...
  , shared_utils.rect.width(displacement), shared_utils.rect.height(displacement) );

use_time_of_first_entry_in_target = params.use_time_of_first_entry_in_target;
allow_out_of_bounds_end_point = params.allow_out_of_bounds_end_point;

%%

t = aligned_outs.params.look_back:aligned_outs.params.look_ahead;
assert( numel(t) == size(aligned_outs.x, 2), 'Time does not match position' );

rts = nan( rows(aligned_outs.labels), 1 );

for i = 1:numel(saccade_out_mask)
  ind = saccade_out_mask(i);
  roi = directional_go_targ_rois(ind, :);
  ib = bfw.bounds.rect( stops(ind, 1), stops(ind, 2), roi );
  
  start_stop_ind = saccade_outs.saccade_start_stop(ind, :);
  aligned_ind = saccade_outs.aligned_to_saccade_ind(ind);
  aligned_x = aligned_outs.x(aligned_ind, start_stop_ind(1):start_stop_ind(2));
  aligned_y = aligned_outs.y(aligned_ind, start_stop_ind(1):start_stop_ind(2));
  all_ib = arrayfun( @(x, y) bfw.bounds.rect(x, y, roi), aligned_x, aligned_y );
  
  if ( ib )
    assert( all_ib(end) );
    
    if ( use_time_of_first_entry_in_target )
      tmp_ind = numel( all_ib );

      while ( tmp_ind > 0 && all_ib(tmp_ind) )
        tmp_ind = tmp_ind - 1;
      end

      offset_ind = tmp_ind + start_stop_ind(1);
      actual_t = t(offset_ind);
    else
      actual_t = t(start_stop_ind(2));
    end
    
    rts(aligned_ind) = actual_t;
  elseif ( allow_out_of_bounds_end_point && any(all_ib) )
    [starts, durs] = shared_utils.logical.find_islands( all_ib );

    if ( ~isempty(starts) )
      if ( use_time_of_first_entry_in_target )
        start_ind = starts(1);
      else
        start_ind = starts(1) + durs(1) - 1;
      end
      
      offset_ind = start_ind + start_stop_ind(1) - 1;
      
      rts(aligned_ind) = t(offset_ind);
    end
  end
end

end