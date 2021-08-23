function [data, labels, I] = maybe_apply_rowop(data, labels, each, op, mask)

%   MAYBE_APPLY_ROWOP -- Conditionally apply row operation to subsets of data.
%
%     [data, labels] = ... maybe_apply_rowop( data, labels, each, op );
%     conditionaly applies function `op` to subsets of `data` identified by
%     unique combinations of labels in `each` categories. `labels` is an
%     fcat object with the same number of rows as `data`. `op` is applied 
%     unless `each` is the empty cell array ({}), in which case `data` and 
%     `labels` are returned as copies.
%
%     [data, labels] = ... maybe_apply_rowop( ..., mask ); operates on rows
%     identified by `mask`, a uint64 index vector. In the case that `op` is
%     not applied, the output `data` and `labels` are still the subsets of
%     the input `data` and `labels` given by `mask`.
%
%     [..., I] = ... maybe_apply_rowop( ... ); additionally returns `I`, a
%     cell array of index vectors with respect to the original input `data`
%     and `labels`. If `op` is applied, then `I` has one element for each
%     row of the output `data`, and each element of `I` is the subset of
%     the input `data` that generated the corresponding row of the output
%     `data`. If `op` is not applied, then `I` has a single element such
%     that `output_data = input_data(I{1}, ...)`.
%
%     Ex // 
%
%     data = fcat.example( 'smalldata' );
%     labels = fcat.example();
%     % Compute a mean separately for each 'dose'
%     [d1, l1, I1] = hwwa.maybe_apply_rowop( data, labels, 'dose', @mean );
%     % Compute a mean separately for each 'dose' and 'roi', only for 
%     % 'scrambled' trials
%     [d2, l2, I2] = hwwa.maybe_apply_rowop( ...
%       data, labels, {'dose', 'roi'}, @mean, find(labels, 'scrambled') );
%     % `each` is empty, so the mean is not computed, and 
%     % `d3` is equivalent to `data(I3{1})`
%     [d3, l3, I3] = hwwa.maybe_apply_rowop( ...
%       data, labels, {}, @mean, find(labels, 'scrambled') );
%
%     See also rowop, fcat/findall

assert_ispair( data, labels );

if ( nargin < 5 )
  mask = rowmask( labels );
end

if ( ~iscell(each) || ~isempty(each) )
  [labels, I] = keepeach( labels', each, mask );
  data = rowop( data, I, op );
else
  labels = labels(mask);
  data = rowref( data, mask );
  I = { mask };
end

assert_ispair( data, labels );

end