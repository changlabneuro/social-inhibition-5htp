function [repeats, non_repeats, curr_inds, prev_inds] = ...
  find_repeats(labels, categories, specificity, mask)

if ( nargin < 4 )
  mask = rowmask( labels );
end

if ( nargin < 3 )
  specificity = {};
end

I = findall_or_one( labels, specificity, mask );

repeats = cell( numel(I), 1 );
non_repeats = cell( size(repeats) );
curr_inds = cell( size(repeats) );
prev_inds = cell( size(repeats) );

for i = 1:numel(I)
  values = categorical( labels, categories, I{i} );
  
  if ( rows(values) <= 1 )
    continue;
  end
  
  current = values(2:end, :);
  prev = values(1:end-1, :);
  
  is_repeat = prev == current;
  is_repeat = all( is_repeat, 2 );
  
  repeat_inds = I{i}(is_repeat) + 1;
  non_repeat_inds = I{i}(~is_repeat) + 1;
  
  repeats{i} = repeat_inds;
  non_repeats{i} = non_repeat_inds;
  curr_inds{i} = I{i}(2:end);
  prev_inds{i} = I{i}(1:end-1);
end

repeats = vertcat( repeats{:} );
non_repeats = vertcat( non_repeats{:} );
prev_inds = vertcat( prev_inds{:} );
curr_inds = vertcat( curr_inds{:} );

end