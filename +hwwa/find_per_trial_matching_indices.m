function find_per_trial_matching_indices(src_labels, dest_labels, src_mask, dest_mask)

if ( nargin < 3 )
  src_mask = rowmask( src_labels );
end
if ( nargin < 4 )
  dest_mask = rowmask( dest_labels );
end

each = { 'unified_filename' };
[src_I, src_selectors] = findall( src_labels, each, src_mask );

for i = 1:numel(src_I)
  src_ind = src_I{i};
  dest_ind = find( dest_labels, src_selectors(:, i), dest_mask );
  
  if ( numel(dest_ind) ~= numel(src_ind) )
    error( 'Source + destination subsets do not match.' );
  end
  
  d = 10;
end

end