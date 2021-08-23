function add_5htp_saline_all_corr(ax, ids, x, y)

sal_ids = ids(hwwa.is_in_scatter_ids(ids, 'saline'));
ser_ids = ids(hwwa.is_in_scatter_ids(ids, '5-htp'));

sal_inds = cat( 1, sal_ids.index );
ser_inds = cat( 1, ser_ids.index );
both_inds = [ sal_inds; ser_inds ];

sal_h = hwwa.add_corr( ax, x(sal_inds), y(sal_inds), 'text_prefix', 'saline: ' );
ser_h = hwwa.add_corr( ax, x(ser_inds), y(ser_inds), 'text_prefix', '5-htp: ' );
all_h = hwwa.add_corr( ax, x(both_inds), y(both_inds), 'text_prefix', 'all: ' );

end