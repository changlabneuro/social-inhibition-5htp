function maybe_save_fig(fig, labels, spec, data_type, subdir, params)

if ( params.do_save )
  shared_utils.plot.fullscreen( fig );
  save_p = hwwa.approach_avoid_data_path( params, 'plots', data_type, subdir );
  dsp3.req_savefig( fig, save_p, labels, spec, params.prefix );
end

end