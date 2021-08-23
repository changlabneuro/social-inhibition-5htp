conf = hwwa.config.load();
conf.PATHS.data_root = '/Volumes/My Passport/NICK/Chang Lab 2016/hww_gng/data';

% files_containing = { 'test_', '1017', '1015', '1012', '1011' };
files_containing = { 'test_' };

hwwa.make_edf_trials( 'files_containing', files_containing ...
  , 'config', conf ...
  , 'event', {'go_target_acquired'} ...
  , 'look_back', [0] ...
  , 'look_ahead', [1e3] ...
  , 'overwrite', true ...
  , 'append', true ...
);

%%

unified_p = hwwa.get_intermediate_dir( 'unified', conf );
edf_trials_p = hwwa.get_intermediate_dir( 'edf_trials', conf );
labels_p = hwwa.get_intermediate_dir( 'labels', conf );
% events_p = hwwa.get_intermediate_dir( 'events', conf );
events_p = hwwa.get_intermediate_dir( 'el_events', conf );

% files_containing = { '1017' };
files_containing = { 'test_', '1017' };

mats = hwwa.require_intermediate_mats( [], edf_trials_p, files_containing );

all_labels = fcat();
all_rois = [];
all_inbounds = [];

sr = 1e3;

for i = 1:numel(mats)
shared_utils.general.progress( i, numel(mats) );

edf_trials_file = shared_utils.io.fload( mats{i} );

un_filename = edf_trials_file.unified_filename;

% 500hz for human files

data_file = shared_utils.io.fload( fullfile(unified_p, un_filename) );
labels_file = shared_utils.io.fload( fullfile(labels_p, un_filename) );
events_file = shared_utils.io.fload( fullfile(events_p, un_filename) );
unified_file = shared_utils.io.fload( fullfile(unified_p, un_filename) );

if ( data_file.opts.STRUCTURE.p_target_left ~= 0 )
  continue;
end

fixation_time = unified_file.opts.STIMULI.setup.go_target.target_duration;

try
  target_event = edf_trials_file.trials('go_target_acquired');
  
  divisor = shared_utils.struct.field_or( target_event, 'sample_rate', 1 );

%   targ_offset_times = events_file.event_times(:, events_file.event_key('go_target_offset'));
  
  targ_acquire_times = events_file.event_times(:, events_file.event_key('go_target_acquired'));
  targ_offset_times = targ_acquire_times + fixation_time;
  
catch err
  continue;
end

amt_ahead = ceil(targ_offset_times - targ_acquire_times) * sr / divisor;
assert( ~any(amt_ahead == 0) );

x = target_event.samples('posX');
y = target_event.samples('posY');

go_targ = unified_file.opts.STIMULI.go_target;

% go_targ.put( 'center-right' );
targ = go_targ.targets{1};

x_off = targ.x_offset;
y_off = targ.y_offset;
offsets = [ x_off, y_off, x_off, y_off ];
right_bounds = targ.bounds + targ.padding + offsets;

p_ib = rowzeros( rows(x) );

for j = 1:rows(x)
  c_amt_ahead = amt_ahead(j);
  
  if ( isnan(c_amt_ahead) ), continue; end
  
  c_amt_ahead = min( c_amt_ahead, size(x, 2) );
  
  subset_x = x(j, 1:c_amt_ahead);
  subset_y = y(j, 1:c_amt_ahead);
  
  is_ib = bfw.bounds.rect( subset_x, subset_y, right_bounds );
  
  p_ib(j) = pnz( is_ib );
end

all_inbounds = [ all_inbounds; p_ib ];

append( all_labels, labels_file.labels );

end

assert_ispair( all_inbounds, all_labels );

hwwa.add_data_set_labels( all_labels );

%%

un_filenames = combs( all_labels, 'unified_filename' );
un_filenames = shared_utils.cell.containing( un_filenames, '1017' );

to_find = cshorzcat( {'correct_true', 'go_trial'}, un_filenames );

is_target_trial = trueat( all_labels, find(all_labels, to_find) );

sum( is_target_trial & all_inbounds > 0.9 ) / sum(is_target_trial) * 100

%%

x = target_event.samples('posX');
y = target_event.samples('posY');

go_targ = data_file.opts.STIMULI.go_target;
go_targ.put('center-right');
right_bounds = go_targ.targets{1}.bounds;
go_targ.put('center-left');
left_bounds = go_targ.targets{1}.bounds;

is_target_trial = find( labels_file.labels, 'correct_false' );

x_incorrect = x(is_target_trial, :);
y_incorrect = y(is_target_trial, :);

all_nans = all( isnan(x_incorrect) & isnan(y_incorrect), 2 );

nan_trial_indices = is_target_trial(all_nans);
