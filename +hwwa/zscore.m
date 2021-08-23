function data = zscore(data, labels, within)

I = findall( labels, within );

colons = repmat( {':'}, 1, ndims(data)-1 );

for i = 1:numel(I)
  subset_data = data(I{i}, colons{:});
  
  means = nanmean( subset_data, 1 );
  devs = nanstd( subset_data, [], 1 );
  
  data(I{i}, colons{:}) = (subset_data - means) ./ devs;
end

end