function result = single_percent_correct(labels, varargin)

n_corr = numel( find(labels, 'correct_true', varargin{:}) );
n_incorr = numel( find(labels, 'correct_false', varargin{:}) );

result = n_corr / (n_corr + n_incorr);

end