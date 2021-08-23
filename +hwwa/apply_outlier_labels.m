function dest_labels = apply_outlier_labels(source_labels, in_bounds, each, dest_labels, src_mask, dest_mask)

if ( nargin < 5 )
  src_mask = rowmask( source_labels );
end

if ( nargin < 6 )
  dest_mask = rowmask( dest_labels );
end

validateattributes( in_bounds, {'logical'}, {'vector', 'numel', rows(source_labels)} ...
  , mfilename, 'in bounds' );

[behav_I, behav_C] = findall( dest_labels, each, dest_mask );
addsetcat( dest_labels, 'is_outlier', 'is_outlier__true' );

for i = 1:numel(behav_I)
  behav_ind = behav_I{i};
  pupil_ind = find( source_labels, behav_C(:, i), src_mask );
  assert( numel(pupil_ind) == numel(behav_ind) );
  
  setcat( dest_labels, 'is_outlier', 'is_outlier__false', behav_ind(in_bounds(pupil_ind)) );
end

end