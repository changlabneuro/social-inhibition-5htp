conf = hwwa.config.load();
conf.PATHS.data_root = '/Volumes/My Passport/NICK/Chang Lab 2016/hww_gng/data';

hwwa.make_first_saccade_info( ...
    'config', conf ...
  , 'debug', false ...
  , 'overwrite', true ...
  , 'pad_beginning', 10 ...
  , 'events', 'go_nogo_cue_onset'...
);

%%

conf = hwwa.config.load();

conf.PATHS.data_root = '/Volumes/My Passport/NICK/Chang Lab 2016/hww_gng/data';

unified_p = hwwa.get_intermediate_dir( 'unified', conf );
labels_p = hwwa.get_intermediate_dir( 'labels', conf );
saccade_p = hwwa.get_intermediate_dir( 'first_saccade/go_nogo_cue_onset', conf );
events_p = hwwa.get_intermediate_dir( 'events', conf );

label_mats = hwwa.require_intermediate_mats( [], labels_p, [] );

summary = labeled();

delays = 0.1:0.01:0.5;
n_delays = 3;
[g_starts, g_stops] = hwwa.bin_delays( delays, n_delays );

g_starts = [ 0.1, 0.2, 0.3, 0.4 ];
g_stops = [0.19, 0.29, 0.39, 0.5 ];

rt_dat = [];
all_saccade_directions = [];
all_is_saccade_before_target_onset = logical( [] );
all_labs = fcat();

for i = 1:numel(label_mats)
  hwwa.progress( i, numel(label_mats) );
  
  labels_file = shared_utils.io.fload( label_mats{i} );
  
  try
    unified_file = shared_utils.io.fload( fullfile(unified_p, labels_file.unified_filename) );
    events_file = shared_utils.io.fload( fullfile(events_p, labels_file.unified_filename) );
    
    saccade_filename = fullfile( saccade_p, labels_file.unified_filename );
    
    if ( shared_utils.io.fexists(saccade_filename) )
      saccade_file = shared_utils.io.fload( saccade_filename );
      has_saccade_file = true;
    else
      warning( 'No saccade file for: "%s".', labels_file.unified_filename );
      has_saccade_file = false;
    end
  catch err
    hwwa.print_fail_warn( labels_file.unified_filename, err.message );
    continue;
  end
  
  go_targ_acq = events_file.event_times(:, events_file.event_key('go_target_acquired'));
  go_cue_on = events_file.event_times(:, events_file.event_key('go_nogo_cue_onset'));
  go_targ_on = events_file.event_times(:, events_file.event_key('go_target_onset'));
  
  targ_offsets = go_targ_acq - go_targ_on;
  
  if ( has_saccade_file )
    saccade_starts = saccade_file.start_times / 1e3;
    directions = saccade_file.directions(:);
  
    is_saccade_before_target_onset = saccade_starts < targ_offsets | isnan(targ_offsets);
  else
    is_saccade_before_target_onset = false( numel(targ_offsets), 1 );
    directions = rownan( numel(targ_offsets) );
  end
  
  rt = [ unified_file.DATA(:).reaction_time ];
  rt = rt(:);
  
  labs = labels_file.labels;
  
  hwwa.add_day_labels( labs );
  hwwa.add_group_delay_labels( labs, g_starts, g_stops );
  hwwa.add_data_set_labels( labs );
  hwwa.add_drug_labels( labs );
  
  addsetcat( labs, 'target-looks', 'no_look_to_target' );
  setcat( labs, 'target-looks', 'looked_to_target', find(~isnan(go_targ_acq)) );
  
  append( all_labs, prune(labs) );
  rt_dat = [ rt_dat; rt ];
  all_saccade_directions = [ all_saccade_directions; directions(:) ];
  all_is_saccade_before_target_onset = [ all_is_saccade_before_target_onset; is_saccade_before_target_onset ];
end

assert_ispair( rt_dat, all_labs );
assert_ispair( all_saccade_directions, all_labs );

is_valid_saccade = ~isnan( all_saccade_directions ) & all_is_saccade_before_target_onset;

plot_p = fullfile( conf.PATHS.data_root, 'plots', 'behavior', datestr(now, 'mmddyy') );

%%  add saccade direction post cue display

addsetcat( all_labs, 'saccade_direction', 'no_saccade' );

is_saccade_left = all_saccade_directions == 0;
is_saccade_right = all_saccade_directions == 1;

setcat( all_labs, 'saccade_direction', 'saccade_left', find(is_saccade_left) );
setcat( all_labs, 'saccade_direction', 'saccade_right', find(is_saccade_right) );

%%  new error types

uselabs = all_labs';

% base_mask = fcat.mask( uselabs ...
%   , @find, hwwa_get_first_days() ...
%   , @findnone, 'no_fixation' ...
% );

base_mask = fcat.mask( uselabs ...
  , @find, hwwa_get_cron_days() ...
  , @findnone, 'no_fixation' ...
);

% [~, errcat] = hwwa_add_new_error_types( uselabs, is_valid_saccade, base_mask );
[~, errcat] = hwwa_add_new_error_types_cron( uselabs, base_mask );

%%  

count_spec = {'cue_delay', 'drug', 'trial_type'};

I = findall( uselabs, count_spec, base_mask );
kinds = combs( uselabs, errcat, base_mask );
n_kinds = numel( kinds );

countlabs = fcat();
counts = [];

stp = 1;

for i = 1:numel(I)
  cs = double( count(uselabs, kinds, I{i}) );
  
  p = cs / numel( I{i} );  
  
  append1( countlabs, uselabs, I{i}, n_kinds );
  setcat( countlabs, errcat, kinds, stp:stp+n_kinds-1 );
  
  counts = [ counts; p(:) ];
  
  stp = stp + n_kinds;
end

%%

do_save = true;

prefix = 'cron';

pltdat = counts;
pltlabs = countlabs';

pl = plotlabeled.make_common();
pl.group_order = { 'is_correct', 'go_never_looked_to_target', 'go_anticipatory_saccade', 'go_broke_cue_fixation' ...
  , 'nogo_looked_to_target', 'nogo_anticipatory_saccade'};

xcats = { 'cue_delay' };
gcats = { errcat };
pcats = { 'drug', 'trial_type' };
pl.panel_order = { '5-htp', 'saline' };

pl.stackedbar( pltdat, pltlabs, xcats, gcats, pcats );

shared_utils.plot.fullscreen( gcf );

if ( do_save )
  dsp3.req_savefig( gcf, plot_p, pltlabs, cshorzcat(pcats), sprintf('%s_summarized_performance', prefix) );
end








