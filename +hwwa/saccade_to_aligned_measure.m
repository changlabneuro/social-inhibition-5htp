function out_measure = saccade_to_aligned_measure(measure, aligned_to_saccade, num_rows, mask)

out_measure = nan( num_rows, 1 );

if ( nargin == 4 )
  % if supplied mask
  masked_to_aligned = aligned_to_saccade(mask);
  out_measure(masked_to_aligned) = measure(mask);
else
  out_measure(aligned_to_saccade) = measure;
end


end