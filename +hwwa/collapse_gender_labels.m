function labs = collapse_gender_labels(labs)

image_categories = incat( labs, 'target_image_category' );

for i = 1:numel(image_categories)
  to_replace = image_categories{i};
  
  male_ind = strfind( lower(to_replace), 'male' );
  female_ind = strfind( lower(to_replace), 'female' );
  
  if ( ~isempty(female_ind) )
    replace_with = to_replace(female_ind + numel('female')+1:end);
  elseif ( ~isempty(male_ind) )
    replace_with = to_replace(male_ind + numel('male')+1:end);
  else
    continue;
  end
  
  replace( labs, to_replace, replace_with );
end

end