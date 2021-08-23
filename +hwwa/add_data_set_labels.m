function labs = add_data_set_labels(labs)

addcat( labs, 'data_set' );

ro1_days = { '050818', '050918', '051018', '051118', '051418' };

tarantino_check_days = { '100118', '100218', '100318', '100418' };

tarantino_check_drug_days = { '101118', '101218', '101518', '101718' ...
  , '102218', '102318' '102418', '102518' };

tarantino_social_days = { '121118', '121218', '121318', '121418', '121718', '121918' };

pilot_social_days = { '011619', '011719', '011819', '012119', '012219' ...
  , '012319', '012519', '012919', '013019', '013119', '020119' ...
  , '020519', '020619', '020719', '020819', '021119', '021319' ...
  , '021419', '021819', '021919', '022119', '022419', '022519' ...
  , '022619', '022719', '022819', '030419', '030519' ...
  , '031819', '031919', '032019', '032119', '032219' ...
  , '041019', '041219', '041619', '041819', '042319', '042519' ...
  , '040819', '041119', '041519', '041719', '042219', '042419' ...
  , '042919', '043019', '050219', '050319' ...
};

add_days( labs, ro1_days, 'ro1' );
add_days( labs, tarantino_check_days, 'tarantino_delay_check' );
add_days( labs, tarantino_check_drug_days, 'tarantino_drug_check' );
add_days( labs, tarantino_social_days, 'tarantino_social' );
add_days( labs, pilot_social_days, 'pilot_social' );

end

function add_days(labs, days, lab)

for i = 1:numel(days)
  setcat( labs, 'data_set', lab, find(labs, days{i}) );
end

end