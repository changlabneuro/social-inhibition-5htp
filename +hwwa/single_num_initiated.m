function result = single_num_initiated(labels, varargin)

result = numel( find(labels, 'initiated_true', varargin{:}) );

end