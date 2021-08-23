function [pupil, I] = median_normalize_pupil(pupil, labels, each, mask)

assert_ispair( pupil, labels );
validateattributes( pupil, {'double'}, {'vector'}, mfilename, 'pupil size' );

if ( nargin < 4 )
  mask = rowmask( labels );  
end

I = findall( labels, each, mask );

for i = 1:numel(I)
  subset = pupil(I{i});
  med = nanmedian( subset );
  normed = subset / med;
  pupil(I{i}) = normed;  
end

end