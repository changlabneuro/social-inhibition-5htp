conf = hwwa_set_local_data_root();

p = fullfile( hwwa.dataroot(conf), 'raw_redux' );

edfs = shared_utils.io.find( p, '.edf' );
edfs = shared_utils.cell.containing( edfs, 'test' );

all_blink_durs = [];

for i = 1:numel(edfs)
  try
    edf_obj = Edf2Mat( edfs{j} );
  catch err
    continue;
  end

  eblink = edf_obj.Events.Eblink;
  blink_start = eblink.start;
  blink_stop = eblink.end;
  
  blink_durs = blink_stop - blink_start;
  
  all_blink_durs = [ all_blink_durs; blink_durs(:) ];
end

%%

hist( all_blink_durs );

med = median( all_blink_durs );

hold on;

shared_utils.plot.add_vertical_lines( gca, med );