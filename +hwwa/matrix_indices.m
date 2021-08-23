function [labels, inds] = matrix_indices(labels, row_spec, column_spec)

%%

column_I = findall( labels, column_spec );
remaining_spec = setdiff( getcats(labels), column_spec );

rows = max( cellfun(@numel, column_I) );
cols = numel( column_I );

mat = nan( rows, cols );
categ = categorical( labels, remaining_spec );
new_categ = categorical();

for i = 1:numel(column_I)
  row_inds = column_I{i};
  
  for j = 1:numel(row_inds)
    row_v = categ(row_inds(j), :);
    
    if ( isempty(new_categ) )
      new_categ(j, :) = row_v;
    else
      ind = all( new_categ == row_v, 2 );
    end
  end
end

end