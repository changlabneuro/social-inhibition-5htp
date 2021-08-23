function save_tabular_data_labels(save_p, data, data_header, labels)

data_header = cellstr( data_header );

validateattributes( data, {'double'}, {'ndims', 2}, mfilename, 'data' );
assert_ispair( data, labels );
assert( numel(data_header) == size(data, 2), 'Header does not match data columns.' );

save_dir = fileparts( save_p );
shared_utils.io.require_dir( save_dir );

[c, label_categories] = categorical( labels );
label_indices = double( c );
label_entries = categories( c );

save( save_p, 'data', 'data_header', 'label_indices', 'label_entries', 'label_categories' );

end