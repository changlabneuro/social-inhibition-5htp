function make_mua_spikes(varargin)

defaults = hwwa.get_common_make_defaults();
defaults.f1 = 700;
defaults.f2 = 1800;
defaults.n = 2;
defaults.std_threshold = 3;

params = hwwa.parsestruct( defaults, varargin );

lfp_p = hwwa.get_intermediate_dir( 'wb' );
output_p = hwwa.get_intermediate_dir( 'mua_spikes' );

f1 = params.f1;
f2 = params.f2;
n = params.n;
thresh = params.std_threshold;

mats = hwwa.require_intermediate_mats( params.files, lfp_p, params.files_containing );

for i = 1:numel(mats)
  hwwa.progress( i, numel(mats), mfilename );
  
  lfp_file = shared_utils.io.fload( mats{i} );
  
  un_filename = lfp_file.unified_filename;
  output_filename = fullfile( output_p, un_filename );
  
  if ( hwwa.conditional_skip_file(output_filename, params.overwrite) )
    continue;
  end
  
  n_chans = size( lfp_file.lfp, 1 );
  
  fs = lfp_file.sample_rate;
  
  first_unit = true;
  
  for j = 1:n_chans
    
    signals = lfp_file.lfp(j, :);
    filtered = hwwa.zpfilter( signals(:), f1, f2, fs, n );
    is_spike = hwwa.get_mua_is_spike( filtered(:), thresh );
    spike_ts = lfp_file.time(is_spike);
    
    channel_name = lfp_file.channel{j};
    unit_id = lfp_file.id{j};
    unit_n = j;
    
    units = struct();
    units.times = spike_ts;
    units.channel = channel_name;
    units.number = unit_n;
    units.id = sprintf( '%s-%d', unit_id, unit_n );

    if ( first_unit )
      all_units = units;
      first_unit = false;
    else
      all_units = [ all_units; units ];
    end
  end
  
  units = struct();
  units.units = all_units;
  units.unified_filename = un_filename;
  
  shared_utils.io.require_dir( output_p );
  
  hwwa.psave( output_filename, units, 'units' );
end

end