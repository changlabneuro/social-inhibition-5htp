function conf = hwwa_set_local_data_root(conf, p)

if ( nargin < 1 || isempty(conf) )
  conf = hwwa.config.load();
end

if ( nargin < 2 || isempty(p) )
  p = '/Volumes/My Passport/NICK/Chang Lab 2016/hww_gng/data';
end

conf.PATHS.data_root = p;

end