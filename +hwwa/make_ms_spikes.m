function make_ms_spikes(varargin)

defaults = hwwa.get_common_make_defaults();
defaults.sample_rate = 40e3;
defaults.mua = false;

params = hwwa.parsestruct( defaults, varargin );

conf = hwwa.config.load();

firings_p = fullfile( conf.PATHS.data_root, 'ms', 'firings' );
meta_p = fullfile( conf.PATHS.data_root, 'ms', 'meta' );
output_p = hwwa.get_intermediate_dir( 'spikes' );

files = params.files;
fc = params.files_containing;

mda_files = hwwa.require_intermediate_mats( files, firings_p, fc, '.mda' );

per_file = containers.Map();

for i = 1:numel(mda_files)
  hwwa.progress( i, numel(mda_files), mfilename );
  
  [~, meta_id] = fileparts( mda_files{i} );
  
  firings = readmda( mda_files{i} );
  meta_file = shared_utils.io.fload( fullfile(meta_p, sprintf('%s.mat', meta_id)) );
  
  un_id = meta_file.unified_id;
  un_filename = meta_file.unified_filename;
  
  output_filename = fullfile( output_p, un_filename );
  
  if ( hwwa.conditional_skip_file(output_filename, params.overwrite) )
    continue;
  end
  
  ms_channel_ids = firings(1, :);
  ms_spike_times = firings(2, :);
  ms_unit_ids = firings(3, :);
  
  if ( params.mua )
    ms_unit_ids = make_mua( ms_unit_ids, ms_channel_ids );
    ms_unit_prefix = 'MS-MUA';
  else
    ms_unit_prefix = 'MS';
  end
  
  unq_channels = unique( ms_channel_ids );
  
  first_unit = true;
  
  for j = 1:numel(unq_channels)
    
    ms_channel_id = unq_channels(j);
    
    channel_ind = ms_channel_ids == ms_channel_id;
    
    channel_name = meta_file.channel_strs{ms_channel_id};
    ms_channel_name = shared_utils.general.channel2str( ms_unit_prefix, ms_channel_id );
    
    unit_ns = unique( ms_unit_ids(channel_ind) );
    
    for k = 1:numel(unit_ns)
      
      unit_n = unit_ns(k);
      unit_ind = channel_ind & ms_unit_ids == unit_n;
      spike_ts = ms_spike_times(unit_ind);
      
      unit = struct();
      unit.times = (spike_ts - 1) ./ params.sample_rate;
      unit.channel = channel_name;
      unit.region = meta_file.region;
      unit.ms_channel_number = ms_channel_id;
      unit.ms_channel = ms_channel_name;
      unit.number = k;
      unit.id = sprintf( '%s-%s-%s-%d', un_id, channel_name, ms_channel_name, k );
      
      if ( first_unit )
        all_units = unit;
        first_unit = false;
      else
        all_units = [ all_units; unit ];
      end
    end
  end
  
  if ( isKey(per_file, output_filename) )
    current = per_file(output_filename);
    current.units = [ current.units; all_units ];
    per_file(output_filename) = current;
  else
    units = struct();
    units.units = all_units;
    units.unified_filename = un_filename;
    per_file(output_filename) = units;
  end
end

K = keys( per_file );

for i = 1:numel(K)
  
  output_filename = K{i};
  units = per_file(output_filename);

  shared_utils.io.require_dir( output_p );
  
  hwwa.psave( output_filename, units, 'units' );
end

end

function unit_ids = make_mua(unit_ids, channel_ids)

unq_channels = unique( channel_ids );

for i = 1:numel(unq_channels)
  ind = channel_ids == unq_channels(i);
  unit_ids(ind) = i;
end

end