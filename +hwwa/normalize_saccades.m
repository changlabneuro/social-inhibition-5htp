function saccades = normalize_saccades(saccades)

for i = 1:size(saccades, 1)
  dir = saccades(i, :);
  saccades(i, :) = dir ./ sqrt( sum(dir .* dir) );
end

end