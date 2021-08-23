function data_labs = merge_unit_labs(data_labs, psth_labs)

import shared_utils.assertions.*;

assert__isa( data_labs, 'fcat', 'data labels' );
assert__isa( psth_labs, 'fcat', 'psth labels' );

n_units = numel( incat(psth_labs, 'id') );
data_sz = size( data_labs, 1 );
psth_sz = size( psth_labs, 1 );

if ( psth_sz ~= data_sz * n_units )
  error( 'PSTH size must be %d times data size.', n_units );
end

repmat( data_labs, n_units );
prune( merge(data_labs, psth_labs) );

end