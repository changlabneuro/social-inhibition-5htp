function [pcorr_dat, pcorr_labs, num_selected] = percent_correct(labels, spec, mask)

if ( nargin < 3 )
  mask = rowmask( labels );
end

[pcorr_labs, I] = keepeach_or_one( labels', spec, mask );
pcorr_dat = rowzeros( numel(I) );
num_selected = rowzeros( numel(I) );

for i = 1:numel(I)
  n_corr = numel( find(labels, 'correct_true', I{i}) );
  n_incorr = numel( find(labels, 'correct_false', I{i}) );
  num_selected(i) = n_corr + n_incorr;

  pcorr_dat(i) = n_corr / (n_corr + n_incorr);
end

end