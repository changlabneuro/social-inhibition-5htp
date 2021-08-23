files = hwwa.approach_avoid_files();

un_files = shared_utils.io.findmat( hwwa.gid('unified') );
un_files = shared_utils.cell.containing( un_files, files );

deg_offs = zeros( numel(un_files), 1 );
px_offs = zeros( size(deg_offs) );

fix_pos = zeros( numel(un_files), 2 );
fix_size = zeros( numel(un_files), 2 );

targ_pos_l = zeros( numel(un_files), 2 );
targ_pos_r = zeros( numel(un_files), 2 );
targ_size = zeros( numel(un_files), 2 );

targ_center = zeros( numel(un_files), 2 );

for i = 1:numel(un_files)
  shared_utils.general.progress( i, numel(un_files) );
  
  un_file = shared_utils.io.fload( un_files{i} );
  off = un_file.config.STIMULI.setup.go_target.displacement;
  w_center = un_file.opts.WINDOW.center;
  w_rect = un_file.opts.WINDOW.rect;
  w_width = shared_utils.rect.width( w_rect );
 
  fix_size(i, :) = un_file.config.STIMULI.setup.fix_square.size;
  fix_pos(i, :) = w_center;
  
  px_off = w_width * 0.25 + abs(off(1));
  px_offs(i) = px_off;
  
  targ_pos_l(i, :) = [w_center(1) - px_off, w_center(2)];
  targ_pos_r(i, :) = [w_center(1) + px_off, w_center(2)];
  targ_size(i, :) = un_file.config.STIMULI.setup.go_target.size;
  
  deg_offs(i) = hwwa.run_px2deg( px_off, [] );
  targ_center(i, :) = shared_utils.rect.center( un_file.opts.STIMULI.go_target.vertices );
end

%%

monitor_info = hwwa.monitor_constants;
monitor_info.monitor_height_cm = 29;

d = hwwa.run_px2deg( 150, [], [], monitor_info );

%%

eccen = hwwa.run_px2deg( px_offs(1), [], [], monitor_info );

%%

monitor_span = hwwa.run_px2deg( w_width, [], [], monitor_info );

%%

use_fix_pos = fix_pos(1, :);
use_fix_sz = fix_size(1, :);
use_targ_pos = targ_center(1, :);
use_targ_sz = targ_size(1, :);
screen_sz = un_file.opts.WINDOW.rect;

ax = gca;
cla( ax );
xlim( ax, [-150, screen_sz(3)+150] );
ylim( ax, [-150, screen_sz(3)+150] );

fix_r = [ use_fix_pos(1)-use_fix_sz(1)*0.5, use_fix_pos(2)-use_fix_sz(2)*0.5, use_fix_sz ];
targ_r = [ use_targ_pos(1)-use_targ_sz(1)*0.5, use_targ_pos(2)-use_targ_sz(2)*0.5, use_targ_sz ];

rectangle( 'position', screen_sz );
rectangle( 'position', fix_r );
rectangle( 'position', targ_r );

axis( ax, 'square' );

