function [num_init, init_labels, each_I] = num_initiated(labels, each, mask)

[init_labels, each_I] = keepeach( labels', each, mask );
num_init = nan( size(each_I) );

for i = 1:numel(each_I)
  num_init(i) = count( labels, 'initiated_true', each_I{i} );
end

end