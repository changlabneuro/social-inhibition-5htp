function pca_outs = saccade_info_pca(saccade_outs, rt, rt_measure_names, perform_each, mask)

validateattributes( rt_measure_names, {'cell'}, {'numel', size(rt, 2)} ...
  , mfilename, 'rt_measure_names' );

[pca_labs, pca_I] = keepeach_or_one( saccade_outs.labels', perform_each, mask );
pca_outs = empty_pca_outs();

for i = 1:numel(pca_I)
  pca_outs(i) = perform_pca( saccade_outs, rt, rt_measure_names, pca_I{i} );
end

end

function s = empty_pca_outs()

fnames = { 'coeff', 'score', 'latent', 'pca_measure_names' };
inputs = cell( 1, numel(fnames)*2 );
inputs(1:2:end) = fnames;
inputs(2:2:end) = {{}};
s = struct( inputs{:} );

end

function factors = make_categorical_factors(num_rows, labels, categories, aligned_to_saccade_ind, mask)

factors = nan( num_rows, numel(categories) );

for i = 1:numel(categories)
  I = findall( labels, categories{i}, mask );
  
  for j = 1:numel(I)
    factors(aligned_to_saccade_ind(I{j}), i) = j-1;
  end
end

end

function pca_outs = perform_pca(saccade_outs, rt, rt_measure_names, mask)

include_categorical_factors = false;
num_rt_measures = size( rt, 2 );

saccade_measure_names = { 'saccade_lengths', 'saccade_velocities' };
% Saccade measures + rt
num_saccade_measures = numel( saccade_measure_names );

num_categorical_factors = ternary( include_categorical_factors, 2, 0 );

pca_measures = nan( rows(rt), num_saccade_measures + num_categorical_factors + num_rt_measures );
pca_measure_names = [ saccade_measure_names, rt_measure_names ];

if ( include_categorical_factors )
  pca_measure_names = [ pca_measure_names, {'drug', 'trial_type'} ];
end

aligned_to_saccade_ind = saccade_outs.aligned_to_saccade_ind;

for i = 1:numel(saccade_measure_names)
  saccade_measure = saccade_outs.(saccade_measure_names{i});
  pca_measures(:, i) = ...
    hwwa.saccade_to_aligned_measure( saccade_measure, aligned_to_saccade_ind, rows(rt), mask );
end

%%

for j = 1:num_rt_measures
  pca_measures(:, j + num_saccade_measures) = rt(:, j);
end

if ( include_categorical_factors )
  pca_measures(:, num_saccade_measures+2:end) = ...
    make_categorical_factors( rows(rt), saccade_outs.labels, {'drug', 'trial_type'} ...
    , saccade_outs.aligned_to_saccade_ind, mask );
end

[coeff, score, latent] = pca( pca_measures );

pca_outs = struct();
pca_outs.coeff = coeff;
pca_outs.score = score;
pca_outs.latent = latent;
pca_outs.pca_measure_names = pca_measure_names;

end