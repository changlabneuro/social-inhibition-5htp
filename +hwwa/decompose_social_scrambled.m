function labs = decompose_social_scrambled(labs)

image_categories = combs( labs, 'target_image_category' );

if ( ~isempty(labs) )
  addsetcat( labs, 'scrambled_type', 'not-scrambled' );
end

for i = 1:numel(image_categories)
  image_category = image_categories{i};
  is_lab = find( labs, image_category );
  
  scr1 = 'scrambled_';
  scr2 = '_scrambled';
  
  scr1_ind = strfind( image_category, scr1 );
  scr2_ind = strfind( image_category, scr2 );
  
  if ( ~isempty(scr1_ind) )
    to_replace = scr1;
    
  elseif ( ~isempty(scr2_ind) )
    to_replace = scr2;
    
  else
    continue;
  end
  
  sans_scrambled = strrep( image_category, to_replace, '' );

  replace( labs, image_category, sans_scrambled );
  setcat( labs, 'scrambled_type', 'scrambled', is_lab );
end

end
