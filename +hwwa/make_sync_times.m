function make_sync_times(varargin)

conf = hwwa.config.load();

defaults = hwwa.get_common_make_defaults();

params = hwwa.parsestruct( defaults, varargin );

conf = params.config;

unified_p = hwwa.get_intermediate_dir( 'unified', conf );
output_p = hwwa.get_intermediate_dir( 'sync', conf );

sync_channel = conf.PLEX.sync_channel;

mats = hwwa.require_intermediate_mats( params.files, unified_p, params.files_containing );

for i = 1:numel(mats)
  hwwa.progress( i, numel(mats), mfilename );
  
  un_file = shared_utils.io.fload( mats{i} );
  
  [pl2_fullfile, pl2_fname] = hwwa.get_pl2_filename( un_file, conf );
  
  if ( ~shared_utils.io.fexists(pl2_fullfile) )
    warning( 'No .pl2 file matching "%s".', pl2_fname );
    continue;
  end
  
  unified_filename = un_file.unified_filename;
  
  output_filename = fullfile( output_p, unified_filename );
  
  if ( hwwa.conditional_skip_file(output_filename, params.overwrite) )
    continue;
  end
  
  sync_vec = PL2Ad( pl2_fullfile, sync_channel );
  
  is_sync_pulse = sync_vec.Values > 4.9e3;
  [pl2_starts, pl2_durs] = shared_utils.logical.find_all_starts( is_sync_pulse(:)' );
  
  pl2_time = (0:sync_vec.FragCounts-1) .* 1/sync_vec.ADFreq;
  
  mat_sync = un_file.opts.PLEX_SYNC;
  mat_sync_times = mat_sync.sync_times(1:mat_sync.sync_iteration-1);
  
  try
    assert( numel(mat_sync_times) == numel(pl2_starts), ['The number of pl2 sync' ...
      , ' times (%d) does not correspond to the number of mat sync times (%d)'] ...
      , numel(mat_sync_times), numel(pl2_starts) );
  catch err
    warning( err.message );
    
    if ( isfield(un_file, 'plex_match_n_sync_times') && un_file.plex_match_n_sync_times )
      
      N = min( numel(mat_sync_times), numel(pl2_starts) );
      
      mat_sync_times = mat_sync_times(1:N);
      pl2_starts = pl2_starts(1:N);
      
    else
      continue;
    end
  end
  
  sync = struct();
  sync.mat = mat_sync_times(:);
  sync.plex = pl2_time(pl2_starts)';
  sync.unified_filename = unified_filename;
  
  shared_utils.io.require_dir( output_p );
  
  hwwa.psave( output_filename, sync, 'sync' ); 
  
end

end