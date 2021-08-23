function [pl2_fullfile, pl2_fname] = get_pl2_filename(un_file, conf)

if ( nargin < 2 )
  conf = hwwa.config.load();
end

data_p = conf.PATHS.data_root;

un_id = un_file.unified_id;
  
pl2_fname = sprintf( '%s.pl2', un_id );

pl2_fullfile = fullfile( data_p, un_file.raw_subdir, pl2_fname );

end
