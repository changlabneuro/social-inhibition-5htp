function [data, labels] = max_normalize(data, labels, spec, varargin)

I = findall( labels, spec, varargin{:} );

for i = 1:numel(I)
  subset = data(I{i}, :, :);
  data(I{i}, :, :) = subset / max( subset(:) );
end

end