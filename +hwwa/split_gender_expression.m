function labs = split_gender_expression(labs)

image_categories = incat( labs, 'target_image_category' );

addcat( labs, 'gender' );

for i = 1:numel(image_categories)
  to_replace = image_categories{i};
  
  male_ind = strfind( lower(to_replace), 'male' );
  female_ind = strfind( lower(to_replace), 'female' );
  
  if ( ~isempty(female_ind) )
    gender = 'female';
    str_ind = female_ind;
    
  elseif ( ~isempty(male_ind) )
    gender = 'male';
    str_ind = male_ind;
    
  else
    continue;
  end
  
  replace_ind = find( labs, to_replace );
  
  replace_with = to_replace(str_ind + numel(gender)+1:end);
  
  replace( labs, to_replace, replace_with );
  setcat( labs, 'gender', gender, replace_ind );
end

end