function rts = saccade_based_reaction_time(aligned_outs, saccade_outs, saccade_out_mask, varargin)

if ( nargin < 3 )
  saccade_out_mask = rowmask( saccade_outs.labels );
end

t = aligned_outs.params.look_back:aligned_outs.params.look_ahead;
assert( numel(t) == size(aligned_outs.x, 2), 'Time does not match position' );

rts = nan( rows(aligned_outs.labels), 1 );

for i = 1:numel(saccade_out_mask)
  ind = saccade_out_mask(i);
  
  start_stop_ind = saccade_outs.saccade_start_stop(ind, :);
  aligned_ind = saccade_outs.aligned_to_saccade_ind(ind);
  
  rts(aligned_ind) = t(start_stop_ind(1));
end

end