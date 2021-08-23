function make_iti_bounds(varargin)

import shared_utils.io.fload;

defaults = hwwa.get_common_make_defaults();
defaults.bin = 1;

params = hwwa.parsestruct( defaults, varargin );

edf_p = hwwa.get_intermediate_dir( 'edf_trials' );
unified_p = hwwa.get_intermediate_dir( 'unified' );
output_p = hwwa.get_intermediate_dir( 'iti_bounds' );

sample_mats = hwwa.require_intermediate_mats( edf_p );

for i = 1:numel(sample_mats)
  hwwa.progress( i, numel(sample_mats), mfilename );
  
  samples = fload( sample_mats{i} );
  
  unified_filename = samples.unified_filename;
  
  unified = fload( fullfile(unified_p, unified_filename) );
  
  output_filename = fullfile( output_p, unified_filename );
  
  if ( hwwa.conditional_skip_file(output_filename, params.overwrite) )
    continue;
  end
  
  iti_x = samples.trials('iti').samples('posX');
  iti_y = samples.trials('iti').samples('posY');
  iti_t = samples.trials('iti').time;
  
  fix_verts = unified.opts.STIMULI.fix_square.vertices;
  targ_padding = unified.opts.STIMULI.fix_square.targets{1}.padding;
  fix_bounds = fix_verts + targ_padding;
  
  ib_x = iti_x >= fix_bounds(1) & iti_x <= fix_bounds(3);
  ib_y = iti_y >= fix_bounds(2) & iti_y <= fix_bounds(4);
  ib = ib_x & ib_y;
  
  if ( params.bin > 1 )
    ib = shared_utils.logical.binned_any( ib, params.bin );
    iti_t = shared_utils.vector.bin( iti_t, params.bin, false );
    iti_t = cellfun( @min, iti_t );
  end
  
  bounds = struct();
  bounds.bounds = ib;
  bounds.time = iti_t;
  bounds.unified_filename = unified_filename;
  bounds.params = params;
  
  shared_utils.io.require_dir( output_p );
  
  save( output_filename, 'bounds' );  
end

end