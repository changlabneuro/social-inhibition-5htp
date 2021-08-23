function make_spikes(varargin)

conf = hwwa.config.load();

defaults = hwwa.get_common_make_defaults();

params = hwwa.parsestruct( defaults, varargin );

unified_p = hwwa.get_intermediate_dir( 'unified' );
output_p = hwwa.get_intermediate_dir( 'spikes' );

mats = hwwa.require_intermediate_mats( params.files, unified_p, params.files_containing );

for i = 1:numel(mats)
  hwwa.progress( i, numel(mats), mfilename );
  
  un_file = shared_utils.io.fload( mats{i} );
  
  [pl2_fullfile, pl2_fname] = hwwa.get_pl2_filename( un_file, conf );
  
  if ( ~shared_utils.io.fexists(pl2_fullfile) )
    warning( 'No .pl2 file matching "%s".', pl2_fname );
    continue;
  end
  
  un_id = un_file.unified_id;
  un_filename = un_file.unified_filename;
  output_filename = fullfile( output_p, un_filename );
  
  if ( hwwa.conditional_skip_file(output_filename, params.overwrite) )
    continue;
  end
  
  pl2_info = PL2GetFileIndex( pl2_fullfile );
  
  spike_chans = pl2_info.SpikeChannels;
  
  active_channels = spike_chans( cellfun(@(x) x.Enabled == 1, spike_chans) );
  
  first_unit = true;
  
  for j = 1:numel(active_channels)
    channel = active_channels{j};
    
    channel_name = channel.Name;
    
    n_units = channel.NumberOfUnits;
    
    for k = 1:n_units
      
      ts = PL2Ts( pl2_fullfile, channel_name, k );
      
      units = struct();
      units.times = ts;
      units.channel = channel_name;
      units.number = k;
      units.id = sprintf( '%s-%s-%d', un_id, channel_name, k );
      
      if ( first_unit )
        all_units = units;
        first_unit = false;
      else
        all_units = [ all_units; units ];
      end
    end
  end
  
  units = struct();
  units.units = all_units;
  units.unified_filename = un_filename;
  
  shared_utils.io.require_dir( output_p );
  
  hwwa.psave( output_filename, units, 'units' );
end

end