function [data, I, ps] = stdthresh(data, labels, within, ndevs)

I = findall( labels, within );

ps = cell( size(I) );

colons = repmat( {':'}, 1, ndims(data)-1 );

for i = 1:numel(I)
  subset_data = data(I{i}, colons{:});

  meaned_data = nanmean( subset_data, 2 );
  
  means = nanmean( meaned_data, 1 );
  devs = nanstd( meaned_data, [], 1 );
  
  amt = devs * ndevs;
  
  lthresh = means - amt;
  uthresh = means + amt;
  
  oob = subset_data < lthresh | subset_data > uthresh;
  
  oob = all( oob, 2 );
  
  subset_data(oob, :) = NaN;
  
  data(I{i}, colons{:}) = subset_data;
  
  ps{i} = oob;
end

end