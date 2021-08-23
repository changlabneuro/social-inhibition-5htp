function out_measures = saccade_to_aligned_measures(saccade_outs, measure_names, num_rows, varargin)

measure_names = cellstr( measure_names );
out_measures = nan( num_rows, numel(measure_names) );

aligned_to_saccade = saccade_outs.aligned_to_saccade_ind;

for i = 1:numel(measure_names)
  measure = saccade_outs.(measure_names{i});
  out_measures(:, i) = ...
    hwwa.saccade_to_aligned_measure( measure, aligned_to_saccade, num_rows, varargin{:} );
end


end