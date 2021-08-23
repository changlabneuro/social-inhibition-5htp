function generate(varargin)

import shared_utils.cell.containing;
import shared_utils.general.json_channels2num;
import shared_utils.general.channel2str;
import shared_utils.cell.ensure_cell;

defaults = hwwa.get_common_make_defaults();
defaults.files_containing = { '.regions' };

params = hwwa.parsestruct( defaults, varargin );

ms_conf = ms_run.config.load();
pipeline_fullfile = fullfile( ms_conf.PATHS.pipelines, 'mountainsort3.mlp' );

conf = hwwa.config.load();
raw_p = fullfile( conf.PATHS.data_root, 'raw_redux' );
ms_p = fullfile( conf.PATHS.data_root, 'ms' );

ms_raw_p = fullfile( ms_p, 'raw' );
ms_script_p = fullfile( ms_p, 'script' );
ms_firings_p = fullfile( ms_p, 'firings' );
ms_meta_p = fullfile( ms_p, 'meta' );

params_fullfile = fullfile( ms_p, 'params', 'params.json' );
ms_runall_fullfile = fullfile( ms_script_p, 'runall.sh' );

files_containing = params.files_containing;

reg_filenames = containing( shared_utils.io.find(raw_p, '.regions'), files_containing );

cellfun( @shared_utils.io.require_dir, {ms_raw_p, ms_script_p, ms_firings_p, ms_meta_p} );

ms_runall_fid = fopen( ms_runall_fullfile, 'wt' );
  
try
  for i = 1:numel(reg_filenames)
    hwwa.progress( i, numel(reg_filenames), mfilename );

    [~, id] = fileparts( reg_filenames{i} );
    pl2_file = sprintf( '%s.pl2', id );
    pl2_fullfile = fullfile( raw_p, pl2_file );

    if ( ~shared_utils.io.fexists(pl2_fullfile) )
      warning( '"%s" does not exist. Skipping ...', pl2_file );
      continue;
    end

    region_map = jsondecode( fileread(reg_filenames{i}) );

    regions = fieldnames( region_map );

    for j = 1:numel(regions)
      region = regions{j};

      chan_ns = json_channels2num( region_map.(region) );
      chans = ensure_cell( channel2str('WB', chan_ns) );

      ms_id = sprintf( '%s_%s', id, region );
      ms_file = sprintf( '%s.mda', ms_id );
      ms_script_file = sprintf( '%s.sh', ms_id );
      ms_meta_file = sprintf( '%s.mat', ms_id );
      
      ms_raw_fullfile = fullfile( ms_raw_p, ms_file );
      ms_firings_fullfile = fullfile( ms_firings_p, ms_file );
      ms_script_fullfile = fullfile( ms_script_p, ms_script_file );
      ms_meta_fullfile = fullfile( ms_meta_p, ms_meta_file );

      ms_run.make_mda_file( pl2_fullfile, chans, ms_raw_fullfile );

      sort_cmd = ms_run.get_sort_command( pipeline_fullfile ...
        , ms_raw_fullfile, ms_firings_fullfile, params_fullfile );

      ms_run_one_fid = fopen( ms_script_fullfile, 'wt' );
      fprintf( ms_run_one_fid, sprintf('%s', sort_cmd) );
      fclose( ms_run_one_fid );
      
      fprintf( ms_runall_fid, sprintf('\n%s', sort_cmd) );
      
      meta_data = struct();
      meta_data.unified_filename = sprintf( '%s.mat', id );
      meta_data.unified_id = id;
      meta_data.region = region;
      meta_data.channel_strs = chans;
      meta_data.channels = chan_ns;
      meta_data.ms_filename = ms_file;
      
      save( ms_meta_fullfile, 'meta_data' );
    end
  end
catch err
  fclose( ms_runall_fid );
  throw( err );
end

end