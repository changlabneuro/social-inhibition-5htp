function tf = is_in_scatter_ids(ids, values)

values = cellstr( values );
tf = arrayfun( @(x) all(ismember(values, x.selectors)), ids );

end