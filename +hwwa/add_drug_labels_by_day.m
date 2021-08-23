function labs = add_drug_labels_by_day(labs)

saline_days = hwwa.get_image_saline_days();
serotonin_days = hwwa.get_image_5htp_days();

add_days( labs, saline_days, 'saline' );
add_days( labs, serotonin_days, '5-htp' );

end

function add_days(labs, days, label_as)

for i = 1:numel(days)
  setcat( labs, 'drug', label_as, find(labs, days{i}) );
end

end