go_targ = hwwa_approach_avoid_go_target_aligned;

%%

[un_I, un_C] = findall( go_targ.labels, 'unified_filename' );

targ_rects = nan( rows(go_targ.labels), 4 );

for i = 1:numel(un_I)
  shared_utils.general.progress( i, numel(un_I) );
  
  un_ind = un_I{i};
  
  trial_dat = hwwa.load1( 'trial_data', un_C{i} );
  unified = hwwa.load1( 'unified', un_C{i} );
  
  screen_center = unified.opts.SCREEN.screens{1}.windows{1}.center;
  screen_size = unified.opts.SCREEN.screens{1}.windows{1}.rect;
  shift = max( unified.config.STIMULI.setup.go_target.displacement );
  targ_size = unified.config.STIMULI.setup.go_target.size;
  
  center_portion_x = (screen_size(3) - screen_size(1)) * 0.25;
  center_left = screen_center(1) - center_portion_x - shift;
  center_right = screen_center(1) + center_portion_x + shift;
  
  left_x0 = center_left - targ_size(1) * 0.5;
  left_x1 = center_left + targ_size(1) * 0.5;
  
  right_x0 = center_right - targ_size(1) * 0.5;
  right_x1 = center_right + targ_size(1) * 0.5;
  
  y0 = screen_center(2) - targ_size(1) * 0.5;
  y1 = screen_center(2) + targ_size(1) * 0.5;
  
  left_p = [ center_left, screen_center(2) ];
  right_p = [ center_right, screen_center(2) ];
  
  left_rect = [ left_x0, y0, left_x1, y1 ];
  right_rect = [ right_x0, y0, right_x1, y1 ];
  
  left_trials = strcmp( {trial_dat.trial_data.target_placement}, 'center-left' );
  assert( numel(left_trials) == numel(un_ind) );
  
  left_inds = un_ind(left_trials);
  for j = 1:numel(left_inds)
    targ_rects(left_inds(j), :) = left_rect;
  end
  right_inds = un_ind(~left_trials);
  for j = 1:numel(right_inds)
    targ_rects(right_inds(j), :) = right_rect;
  end
end

%%

x = go_targ.x;
y = go_targ.y;

targ_dists = nan( size(x) );

for i = 1:numel(un_I)
  shared_utils.general.progress( i, numel(un_I) );
  
  un_ind = un_I{i};
  r = targ_rects(un_ind, :);
  
  for j = 1:numel(un_ind)
    p = [ mean(r(j, [1, 3])), mean(r(j, [2, 4])) ];
    
    for k = 1:size(x, 2)
      px = x(un_ind(j), k);
      py = y(un_ind(j), k);
      targ_dists(un_ind(j), k) = bfw.distance( px, py, p(1), p(2) );
    end
  end
end

%%

to_save = go_targ;
to_save.x = nanmean( to_save.x, 2 );
to_save.y = nanmean( to_save.y, 2 );
to_save.pupil = nanmean( to_save.pupil, 2 );
to_save.mean_targ_dist = nanmean( targ_dists, 2 );

%%

save('/Volumes/My Passport/NICK/Chang Lab 2016/hww_gng/data/intermediates/processed/behavior/go_target_aligned.mat','to_save')