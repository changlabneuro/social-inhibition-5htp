function [normed, norm_labels] = saline_normalize(data, labels, spec, varargin)

non_drug_spec = cssetdiff( spec, 'drug' );
non_drug_I = findall_or_one( labels, non_drug_spec, varargin{:} );

normed = [];
norm_labels = fcat();

for i = 1:numel(non_drug_I)
  sal_ind = find( labels, 'saline', non_drug_I{i} );
  serotonin_ind = find( labels, '5-htp', non_drug_I{i} );
  
  denom = nanmean( data(sal_ind, :) );
  num = data(serotonin_ind, :);
  
  tmp_normed = num ./ denom;
  
  normed = [ normed; tmp_normed ];
  append( norm_labels, labels, serotonin_ind );
end

if ( ~isempty(norm_labels) )
  setcat( norm_labels, 'drug', '5-htp/saline' );
end

assert_ispair( normed, norm_labels );

end