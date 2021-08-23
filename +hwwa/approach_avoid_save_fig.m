function save_p = approach_avoid_save_fig(fig, axs, labels, cats, subdir, params)

shared_utils.plot.fullscreen( fig );
try
  shared_utils.plot.match_ylims( axs );
catch err
  warning( err.message );
end
save_p = hwwa.approach_avoid_data_path( params, 'plots', subdir );
dsp3.req_savefig( fig, save_p, labels, cats, params.prefix );

end