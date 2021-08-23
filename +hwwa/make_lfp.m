function make_lfp(varargin)

conf = hwwa.config.load();

defaults = hwwa.get_common_make_defaults();
defaults.kind = 'fp';

params = hwwa.parsestruct( defaults, varargin );

kind = params.kind;

switch ( kind )
  case 'fp'
    lfp_dir = 'lfp';
    channel_prefix = 'FP';
  case 'wb'
    lfp_dir = 'wb';
    channel_prefix = 'WB';
  otherwise
    error( 'Unrecognized kind "%s".', kind );
end

unified_p = hwwa.get_intermediate_dir( 'unified' );
output_p = hwwa.get_intermediate_dir( lfp_dir );

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
  analog_chans = pl2_info.AnalogChannels;
  
  is_active_spk = cellfun( @(x) x.Enabled == 1, spike_chans );
  active_spk = spike_chans(is_active_spk);
  active_ns = cellfun( @(x) x.Channel, active_spk );
  
  is_active_func = @(x) any(active_ns == x.Channel) && strcmp(x.SourceName, channel_prefix);
  
  is_active_fp = cellfun( is_active_func, analog_chans );
  
  active_channels = analog_chans(is_active_fp);
  
  lfp = struct();
  lfp.unified_filename = un_filename;
  lfp.channel = cell( numel(active_channels), 1 );
  lfp.id = cell( size(lfp.channel) );
  
  for j = 1:numel(active_channels)
    fprintf( '\n\t %d of %d', j, numel(active_channels) );
    
    channel = active_channels{j};
    
    channel_name = channel.Name;
    
    values = PL2Ad( pl2_fullfile, channel_name );
    
    if ( j == 1 )
      all_values = nan( numel(active_channels), numel(values.Values) );
      sample_rate = values.ADFreq;
      t = (0:size(all_values, 2)-1) .* 1/sample_rate;
    end
    
    all_values(j, :) = values.Values;    
    
    lfp.channel{j} = channel_name;
    lfp.id{j} = sprintf( '%s-%s', un_id, channel_name );
  end
  
  lfp.lfp = all_values;
  lfp.time = t;
  lfp.sample_rate = sample_rate;
  
  shared_utils.io.require_dir( output_p );
  
  hwwa.psave( output_filename, lfp, 'lfp', '-v7.3' );
end

end