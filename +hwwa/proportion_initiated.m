function result = proportion_initiated(labels, numerator_mask, denom_mask_func)

% n_init = numel( find(labels, 'initiated_true', varargin{:}) );
% n_non_init = numel( find(labels, 'initiated_false', varargin{:}) );

result = n_init / (n_init + n_non_init);

end