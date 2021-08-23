function saccade_info_pca_plot(pca_outs, labels, groups)

if ( nargin < 3 )
  groups = {};
end

for i = 1:numel(pca_outs)
  ax = gca;
  
  biplot_vectors( ax, pca_outs(i) );
  scatter_groups( ax, pca_outs(i), labels, groups );
end

end

function [data, labels] = run_level_mean(data, labels, varargin)

[labels, I] = keepeach( labels', {'date', 'trial_type', 'drug'}, varargin{:} );
data = bfw.row_nanmean( data, I );

end

function [colsign, max_coeff_len] = coeff_info(coeffs)

% help biplot

[p, d] = size( coeffs );
[~,maxind] = max(abs(coeffs),[],1);
colsign = sign(coeffs(maxind + (0:p:(d-1)*p)));
max_coeff_len = sqrt(max(sum(coeffs.^2, 2)));

end

function scatter_groups(ax, pca_outs, labels, groups)

scores = pca_outs.score;
[colsign, max_coeff_len] = coeff_info( pca_outs.coeff );

[scores, labels] = run_level_mean( scores, labels );

[group_ind, group_combs] = findall( labels, groups );
group_vec = nan( rows(labels), 1 );

for j = 1:numel(group_ind)
  group_vec(group_ind{j}) = j;
end

scores = bsxfun(@times, max_coeff_len.*(scores ./ max(abs(scores(:)))), colsign);

hold( ax, 'on' );
h = gscatter( scores(:, 1), scores(:, 2), group_vec );
legend( h, strrep(fcat.strjoin(group_combs, ' | '), '_', ' ') );

end

function biplot_vectors(ax, pca_outs)

coeff = pca_outs.coeff;
names = pca_outs.pca_measure_names;
scores = pca_outs.score;

hold( ax, 'off' );
biplot( coeff(:, 1:2) ...
  , 'varlabels', strrep(names, '_', ' ') ...
);

end