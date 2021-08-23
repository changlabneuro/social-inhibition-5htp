function ind = find_many(labels, cs, varargin)

ind = [];

for i = 1:size(cs, 1)
  ind = union( ind, find(labels, cs(i, :), varargin{:}) );
end

end