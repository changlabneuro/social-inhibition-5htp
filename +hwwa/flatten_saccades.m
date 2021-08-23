function [saccades, flat_labs, ind, rev_ind] = flatten_saccades(saccades, labels, max_num, mask)

assert_ispair( saccades, labels );

if ( nargin < 4 )
  mask = rowmask( labels );
end

if ( nargin < 3 || isempty(max_num) )
  max_num = inf;
end

trial_counts = cellfun( @rows, saccades );
flat_labs = fcat.like( labels );

subset_counts = trial_counts(mask);
ind = nan( numel(subset_counts), 1 );
rev_ind = zeros( rows(labels), 1 );
num_ind = 0;

for i = 1:numel(subset_counts)
  num_to_append = min( max_num, subset_counts(i) );
  
  for j = 1:num_to_append
    append( flat_labs, labels, mask(i) );
    num_ind = num_ind + 1;
    ind(num_ind) = mask(i);
  end
end

ind = ind(1:num_ind);

saccades = cellfun( @(x) conditional(@() rows(x) < max_num ...
  , x(1:rows(x), :) ...
  , @() x(1:max_num, :)), saccades(mask), 'un', 0 );

saccades = vertcat( saccades{:} );

end