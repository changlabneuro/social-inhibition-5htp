function make_raw_power(varargin)

defaults = hwwa.get_common_make_defaults();
defaults.step_size = 0.05;
defaults.f1 = 2.5;
defaults.f2 = 250;
defaults.n = 2;
defaults.filter = true;
defaults.subtract_reference = true;

params = hwwa.parsestruct( defaults, varargin );

aligned_p = hwwa.get_intermediate_dir( 'aligned_lfp' );
output_p = hwwa.get_intermediate_dir( 'raw_power' );

mats = hwwa.require_intermediate_mats( params.files, aligned_p, params.files_containing );

step_size = params.step_size;

parfor i = 1:numel(mats)
  hwwa.progress( i, numel(mats), mfilename );
  
  lfp = shared_utils.io.fload( mats{i} );
  
  output_filename = fullfile( output_p, lfp.unified_filename );
  
  if ( hwwa.conditional_skip_file(output_filename, params.overwrite) )
    continue;
  end
  
  w_size = lfp.params.window_size;
  
  evts = keys( lfp.psth );
  
  power = containers.Map();
  
  for j = 1:numel(evts)
    fprintf( '\n\t %d of %d', j, numel(evts) );
    
    psth = lfp.psth( evts{j} );
    
    addcat( psth.labels, 'unified_filename' );
    setcat( psth.labels, 'unified_filename', lfp.unified_filename );
    hwwa.add_region_labels( psth.labels, 'FP' );
    
    if ( params.subtract_reference )
      [data, labs] = ref_subtract( psth.data, psth.labels );
      psth.data = data;
      psth.labels = labs;
    end
    
    I = findall( psth.labels, 'id' );
    
    all_chans = [];
    
    for k = 1:numel(I)
      fprintf( '\n\t\t %d of %d', k, numel(I) );
      [one_chan, f, t] = per_channel( psth, params, w_size, step_size, I{k} );
      
      all_chans = [ all_chans; one_chan ];
    end
    
    power( evts{j} ) = struct( 'data', all_chans, 'labels', psth.labels ...
      , 'time', t, 'frequencies', f, 'window_size', w_size ...
      , 'step_size', step_size );
  end
  
  lfp = rmfield( lfp, 'psth' );
  lfp.measure = power;
  lfp.params = params;
  
  shared_utils.io.require_dir( output_p );
  
  fprintf( '\n Saving ... ' );
  hwwa.psave( output_filename, lfp, 'power', '-v7.3' );
  fprintf( 'Done.' );
end

end

function [data, labs] = ref_subtract( data, labs )

ref_ind = find( labs, 'ref' );

if ( isempty(ref_ind) )
  return;
end

ref_data = data(ref_ind, :);

I = findall( labs, 'channel' );

for i = 1:numel(I)
  subset_data = data(I{i}, :);
  data(I{i}, :) = subset_data - ref_data;
end

end

function [all_c, f, t] = per_channel(psth, params, window_size, step_size, ind)

colons = repmat( {':'}, 1, ndims(psth.data)-1 );

data = psth.data(ind, colons{:});

data = data * 1e4;

sr = psth.sample_rate;

if ( params.filter )
  f1 = params.f1;
  f2 = params.f2;
  n = params.n;      
  data = hwwa.zpfilter( data, f1, f2, sr, n );
end

win_samples = window_size * sr;
stp_samples = step_size * sr;

data = shared_utils.array.bin3d( data, win_samples, stp_samples );

chronux_params = struct();
chronux_params.Fs = sr;
chronux_params.tapers = [ 1.5, 2 ];

t = psth.time(1:stp_samples:end-win_samples+1);

for i = 1:size(data, 3)
  
  one_window = data(:, :, i);
  
  [pxx, f] = periodogram( one_window', [], 1:200, sr );
  
  pxx = pxx';
  
%   for k = 1:size(one_window, 1)
% %     [pow, f] = mtspectrumc( one_window(k, :), chronux_params );
%     [pow, f] = periodogram( one_window(k, :)', [], 1:200, sr );
%     
%     if ( k == 1 )
%       pxx = nan( size(one_window, 1), numel(pow) );
%     end
%     
%     pxx(k, :) = pow;
%   end

  if ( i == 1 )
    all_c = nan( [size(pxx), size(data, 3)] );
  end

  all_c(:, :, i) = pxx;
end

end