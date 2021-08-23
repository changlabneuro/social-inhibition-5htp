function samples_file = edf_samples(files)

hwwa.validatefiles( files, 'edf' );

edf_file = shared_utils.general.get( files, 'edf' );

samples_file = struct();
samples_file.unified_filename = edf_file.unified_filename;

samples_file.t = get_time( edf_file.Samples.time );
samples_file.x = single( edf_file.Samples.posX );
samples_file.y = single( edf_file.Samples.posY );
samples_file.pupil = single( edf_file.Samples.pupilSize );

end

function t = get_time(edf_t)

uint32_max = double( intmax('uint32') );

if ( numel(edf_t) < uint32_max && max(edf_t) < uint32_max && min(edf_t) >= 0 )
  t = uint32( edf_t );
else
  t = edf_t;
end

end