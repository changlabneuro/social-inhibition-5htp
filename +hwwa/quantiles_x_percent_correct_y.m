function [pcorr, pcorr_labels, match_inds] = ...
  quantiles_x_percent_correct_y(data, labels, each, of, num_quants, mask)

quants = dsp3.quantiles_each( data, labels', num_quants, each, of, mask );

full_spec = union( each, of );
each_I = findall( labels, full_spec, intersect(mask, find(~isnan(quants))) );

match_inds = cell( size(each_I) );
pcorr = cell( size(each_I) );
pcorr_labels = cell( size(each_I) );

parfor i = 1:numel(each_I)
  sub_quants = quants(each_I{i});  
  vs = unique( sub_quants );
  
  tmp_pcorr = [];
  tmp_pcorr_labels = fcat();
  tmp_match_inds = {};
  
  for j = 1:numel(vs)
    sub_mask = each_I{i}(sub_quants == vs(j));
    
    [tmp_pcorr(end+1, 1), tmp_labs] = ...
      hwwa.percent_correct( labels', full_spec, sub_mask );
    append( tmp_pcorr_labels, tmp_labs );
    
    tmp_match_inds{end+1, 1} = sub_mask;
  end
  
  quant_labs = arrayfun( @(x) sprintf('x-quantile__%d', x), vs, 'un', 0 );
  addsetcat( tmp_pcorr_labels, 'x-quantile', quant_labs );
  
  pcorr{i} = tmp_pcorr;
  pcorr_labels{i} = tmp_pcorr_labels;
  match_inds{i} = tmp_match_inds;
end

pcorr = vertcat( pcorr{:} );
pcorr_labels = vertcat( fcat, pcorr_labels{:} );
match_inds = vertcat( match_inds{:} );

end