function labs = add_drug_labels(labs)

import shared_utils.assertions.*;

assert__isa( labs, 'fcat', 'labels' );

saline_days = { '0215', '0221', '0228', '0302', '0303', '0305' ...
  , '0508', '0510', '1011', '1017', '1022', '1024' };

serotonin_days = { '0214', '0216', '0220', '0222', '0301', '0304' ...
  , '0509', '0511', '0514', '1012', '1015', '1025', '1023' };

none_days = { '0412', '0410' };

add_labs( labs, saline_days, 'saline' );
add_labs( labs, serotonin_days, '5-htp' );
add_labs( labs, none_days, 'no drug' );

end

function add_labs(labs, days, drug)

for i = 1:numel(days)
  ind = find( labs, [ days{i}, '.mat' ] );
  if ( isempty(ind) ), continue; end
  setcat( labs, 'drug', drug, ind );
end

end