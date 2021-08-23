function labs = fix_image_category_labels(labs)

category = 'target_image_category';

if ( ~hascat(labs, category) )
  return
end

mispellings = { 'appetative', 'scambled', 'lip_smack' };
correct_spellings = { 'appetitive', 'scrambled', 'appetitive' };

assert( numel(mispellings) == numel(correct_spellings) );

labels = incat( labs, category );

for i = 1:numel(mispellings)
  
  has_mispelling = cellfun( @(x) ~isempty(strfind(x, mispellings{i})), labels );
  
  for j = 1:numel(has_mispelling)
    if ( ~has_mispelling(j) )
      continue
    end
    
    current_label = labels{j};
    replace_with = strrep( current_label, mispellings{i}, correct_spellings{i} );
    
    replace( labs, current_label, replace_with );
    
    labels{j} = replace_with;
  end
end

end