function [normed, norm_labels] = scrambled_normalize(data, labels, spec, varargin)

non_scrambled_spec = cssetdiff( spec, 'scrambled_type' );
non_drug_I = findall_or_one( labels, non_scrambled_spec, varargin{:} );

normed = [];
norm_labels = fcat();

for i = 1:numel(non_drug_I)
  scrambled_ind = find( labels, 'scrambled', non_drug_I{i} );
  social_ind = find( labels, 'not-scrambled', non_drug_I{i} );
  
  denom = nanmean( data(scrambled_ind, :) );
  num = data(social_ind, :);
  
  tmp_normed = num ./ denom;
  
  normed = [ normed; tmp_normed ];
  append( norm_labels, labels, social_ind );
end

if ( ~isempty(norm_labels) )
  setcat( norm_labels, 'scrambled_type', 'scrambled/not-scrambled' );
end

assert_ispair( normed, norm_labels );

end