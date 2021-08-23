function labs = tidy_up_labels(labs, monk_label)

if ( nargin < 2 )
  monk_label = '';
end

hwwa.add_day_labels( labs );
hwwa.add_data_set_labels( labs );
hwwa.add_drug_labels_by_day( labs );
hwwa.fix_image_category_labels( labs );
hwwa.split_gender_expression( labs );
hwwa.decompose_social_scrambled( labs );

addcat( labs, 'monkey' );

if ( ~isempty(labs) && ~isempty(monk_label) )
  setcat( labs, 'monkey', lower(monk_label) );
end

end