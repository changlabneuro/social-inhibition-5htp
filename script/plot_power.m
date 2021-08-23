conf = hwwa.config.load();

label_p = hwwa.get_intermediate_dir( 'labels' );
power_p = hwwa.get_intermediate_dir( 'raw_power' );

power_mats = hwwa.require_intermediate_mats( power_p );

evt = 'go_target_onset';
% evt = 'go_target_acquired';
% evt = 'go_nogo_cue_onset';
% evt = 'reward_onset';
% evt = 'go_target_offset';

do_norm = false;
baseline_evt = 'go_nogo_cue_onset';
pre_baseline = -0.3;
post_baseline = 0;
% norm_func = @minus;
norm_func = @rdivide;

do_z = true;
z_within = { 'id' };

psth_labs = fcat();
psth_data = [];

for i = 1:numel(power_mats)
  hwwa.progress( i, numel(power_mats) );
  
  power_file = shared_utils.io.fload( power_mats{i} );
  un_filename = power_file.unified_filename;
  labs_file = shared_utils.io.fload( fullfile(label_p, un_filename) );
  
  measure = power_file.measure(evt);
  psth_d = measure.data;
  
  hwwa.merge_unit_labs( labs_file.labels, measure.labels );
  hwwa.add_broke_cue_labels( labs_file.labels );
  
  if ( do_z )
    psth_d = hwwa.zscore( psth_d, labs_file.labels, z_within );
  end
  
  if ( do_norm )
    baseline_p = power_file.measure(baseline_evt);
    baseline_t = baseline_p.time;
    baseline_ind = baseline_t >= pre_baseline & baseline_t <= post_baseline;
    baseline_d = nanmean( baseline_p.data(:, :, baseline_ind), 3 );
    
    for j = 1:size(psth_d, 3)
      psth_d(:, :, j) = norm_func( psth_d(:, :, j), baseline_d );
    end
  end
  
  append( psth_labs, labs_file.labels );
  psth_data = [ psth_data; psth_d ];
  
  t = measure.time;
end

prune( psth_labs );

save_p = fullfile( conf.PATHS.data_root, 'plots', 'raw_power' ...
  , datestr(now, 'mmddyy'), evt );

%%

cont = SignalContainer( psth_data, SparseLabels.from_fcat(psth_labs) );

cont.frequencies = measure.frequencies;
cont.start = t(1) * 1e3;
cont.stop = t(end) * 1e3;
cont.window_size = measure.window_size * 1e3;
cont.step_size = measure.step_size * 1e3;

%%  correct - incorrect

do_save = false;

figure(1);

specificity = { 'trial_type', 'trial_outcome', 'drug', 'region', 'type_changed' };

meaned = each1d( cont, specificity, @rowops.nanmean );
meaned = rm( meaned, {'no_choice', 'no drug'} );
meaned = collapse( meaned, 'error' );
meaned = only( meaned, 'dlpfc' );

% meaned = only( meaned, 'type_changed_false' );

sub1 = meaned({'go_trial', 'go_choice'}) - meaned({'go_trial', 'nogo_choice'});
sub2 = meaned({'nogo_trial', 'nogo_choice'}) - meaned({'nogo_trial', 'go_choice'});
meaned = [ sub1; sub2 ];

meaned = add_field( meaned, 'epochs', strrep(evt, '_', ' ') );

meaned.spectrogram( {'trial_type', 'trial_outcome', 'drug', 'type_changed'} ...
  , 'frequencies', [0, 100] ...
  , 'time', [-500, 500] ...
  , 'shape', [4, 2] ...
  );

if ( do_save )
  shared_utils.io.require_dir( save_p );

  fname = strjoin( flat_uniques(meaned, {'trial_type', 'trial_outcome'}), '_' );
  full_fname = fullfile( save_p, fname );
  shared_utils.plot.save_fig( gcf, full_fname, {'epsc', 'png', 'fig'}, true );
end

%%  changed minus no-change

do_save = false;

figure(1);

specificity = { 'correct', 'drug', 'region', 'type_changed' };

meaned = each1d( cont, specificity, @rowops.nanmean );
meaned = rm( meaned, {'no_choice', 'no drug'} );
meaned = collapse( meaned, 'error' );
meaned = only( meaned, 'dlpfc' );

% meaned = only( meaned, 'type_changed_false' );

% meaned = only( meaned, 'correct_true' );

sub1 = meaned({'type_changed_true'}) - meaned({'type_changed_false'});

meaned = sub1;

meaned = add_field( meaned, 'epochs', strrep(evt, '_', ' ') );

meaned.spectrogram( specificity ...
  , 'frequencies', [0, 100] ...
  , 'time', [-500, 500] ...
  , 'shape', [2, 2] ...
  );

if ( do_save )
  shared_utils.io.require_dir( save_p );

  fname = strjoin( flat_uniques(meaned, {'trial_type', 'trial_outcome'}), '_' );
  full_fname = fullfile( save_p, fname );
  shared_utils.plot.save_fig( gcf, full_fname, {'epsc', 'png', 'fig'}, true );
end

%%  all 4 trial outcomes

do_save = true;

figure(1);

subset_cont = only( cont, 'initiated_true' );

% subset_cont = rm( subset_cont, '0514.mat' );

specificity = { 'trial_type', 'correct', 'drug', 'region' };

meaned = each1d( subset_cont, specificity, @rowops.nanmean );

meaned = rm( meaned, {'no_choice', 'no drug'} );
meaned = collapse( meaned, 'error' );
meaned = only( meaned, 'dlpfc' );

% meaned = only( meaned, 'saline' );

sub1 = meaned({'go_trial', 'go_choice'}) - meaned({'go_trial', 'nogo_choice'});
sub2 = meaned({'nogo_trial', 'nogo_choice'}) - meaned({'nogo_trial', 'go_choice'});
meaned = [ sub1; sub2 ];

% meaned = meaned({'go_trial'}) - meaned({'nogo_trial'});
% shp = [ 1, 2 ];

meaned = add_field( meaned, 'epochs', strrep(evt, '_', ' ') );

meaned.spectrogram( specificity ...
  , 'frequencies', [0, 65] ...
  , 'time', [-250, 500] ...
  , 'shape', [2, 2] ...
  );

if ( do_save )
  shared_utils.io.require_dir( save_p );

  fname = strjoin( flat_uniques(meaned, {'trial_type', 'trial_outcome'}), '_' );
  full_fname = fullfile( save_p, fname );
  shared_utils.plot.save_fig( gcf, full_fname, {'epsc', 'png', 'fig'}, true );
end

%%  lines

specificity = { 'trial_type', 'correct', 'drug', 'region', 'channel', 'date' };

meaned = each1d( subset_cont, specificity, @rowops.nanmean );

meaned = rm( meaned, {'no_choice', 'no drug'} );
meaned = collapse( meaned, 'error' );
meaned = only( meaned, 'dlpfc' );

sub1 = meaned({'go_trial', 'go_choice'}) - meaned({'go_trial', 'nogo_choice'});
sub2 = meaned({'nogo_trial', 'nogo_choice'}) - meaned({'nogo_trial', 'go_choice'});
meaned = [ sub1; sub2 ];

time_meaned = time_mean( meaned, [-100, 300] );

%%

do_save = true;

pl = ContainerPlotter();
pl.x = time_meaned.frequencies;
pl.main_line_width = 1;
pl.add_ribbon = true;
pl.params.error_function = @plotlabeled.nansem;
pl.compare_series = true;
pl.p_correct_type = 'fdr';
pl.params.smooth_function = @(x) smooth(x, 5);
pl.params.add_smoothing = true;
pl.y_lim = [-0.4, 0.4];

lines_are = { 'drug' };
panels_are = { 'trial_type' };

axs = pl.plot( time_meaned, lines_are, panels_are );

xlim( axs, [0, 65] );

if ( do_save )
  full_save_p = fullfile( save_p, 'lines' );
  
  shared_utils.io.require_dir( full_save_p );

  fname = strjoin( flat_uniques(time_meaned, {'trial_type', 'trial_outcome'}), '_' );
  full_fname = fullfile( full_save_p, fname );
  shared_utils.plot.save_fig( gcf, full_fname, {'epsc', 'png', 'fig'}, true );
end

%%  prev trial correct

do_save = false;

figure(1);

subset_cont = only( cont, 'prevcorrect_true' );

specificity = { 'correct', 'trial_type', 'type_changed', 'drug', 'region' };

meaned = each1d( subset_cont, specificity, @rowops.nanmean );
ns = each1d( subset_cont, specificity, @(x) size(x, 1) );

meaned = rm( meaned, {'no_choice', 'no drug'} );
meaned = collapse( meaned, 'error' );
meaned = only( meaned, 'dlpfc' );

ns = rm( ns, {'no_choice', 'no drug'} );
ns = collapse( ns, 'error' );
ns = only( ns, 'dlpfc' );
ns = only( ns, {'correct_true', 'saline'} );

meaned = add_field( meaned, 'epochs', strrep(evt, '_', ' ') );

meaned = only( meaned, {'correct_true', 'saline'} );

% meaned = meaned({'type_changed_true'}) - meaned({'type_changed_false'});
% meaned = meaned({'nogo_trial'}) - meaned({'go_trial'});
%%
meaned.spectrogram( specificity ...
  , 'frequencies', [0, 100] ...
  , 'time', [-300, 500] ...
  , 'shape', [2, 2] ...
  );

if ( do_save )
  shared_utils.io.require_dir( save_p );

  fname = strjoin( flat_uniques(meaned, {'trial_type', 'trial_outcome'}), '_' );
  full_fname = fullfile( save_p, fname );
  shared_utils.plot.save_fig( gcf, full_fname, {'epsc', 'png', 'fig'}, true );
end


%%  go - nogo

meaned = each1d( cont, {'trial_type', 'trial_outcome', 'drug', 'region'}, @rowops.nanmean );
meaned = rm( meaned, {'no_choice', 'no drug'} );
meaned = collapse( meaned, 'error' );
meaned = only( meaned, 'dlpfc' );

sub1 = meaned({'go_trial', 'go_choice'}) - meaned({'nogo_trial', 'nogo_choice'});
sub2 = meaned({'nogo_trial', 'go_choice'}) - meaned({'go_trial', 'nogo_choice'});

sub1( 'trial_type' ) = 'correct_go minus correct nogo';
sub2( 'trial_type' ) = 'incorrect_go minus incorrect nogo';

meaned = [ sub1; sub2 ];

meaned = add_field( meaned, 'epochs', strrep(evt, '_', ' ') );

meaned.spectrogram( {'trial_type', 'trial_outcome', 'drug'} ...
  , 'frequencies', [0, 100] ...
  , 'shape', [2, 2] ...
  );


