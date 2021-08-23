function add_normalized_corr(ax, ids, x, y, varargin)

ser_ids = ids(hwwa.is_in_scatter_ids(ids, '5-htp/saline'));
ser_inds = cat( 1, ser_ids.index );
ser_h = hwwa.add_corr( ax, x(ser_inds), y(ser_inds), 'text_prefix', '5-htp/saline: ', varargin{:} );

end