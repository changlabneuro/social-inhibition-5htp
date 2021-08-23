function outs = hwwa_check_problematic_go_trials(conf, varargin)

defaults = hwwa.get_common_make_defaults();

params = hwwa.parsestruct( defaults, varargin );

if ( nargin < 1 || isempty(conf) )
  conf = hwwa.config.load();
end

unified_p = hwwa.get_intermediate_dir( 'unified', conf );
edf_trials_p = hwwa.get_intermediate_dir( 'edf_trials', conf );
events_p = hwwa.get_intermediate_dir( 'edf_events', conf );
labels_p = hwwa.get_intermediate_dir( 'labels', conf );

% files_containing = { '1017' };
files_containing = params.files_containing;

mats = hwwa.require_intermediate_mats( [], edf_trials_p, files_containing );

all_labels = fcat();
all_inbounds = [];
all_nan = [];
all_dists = [];
all_x_dists = [];
all_devs = [];

for i = 1:numel(mats)
shared_utils.general.progress( i, numel(mats) );

edf_trials_file = shared_utils.io.fload( mats{i} );

un_filename = edf_trials_file.unified_filename;

labels_file = shared_utils.io.fload( fullfile(labels_p, un_filename) );
events_file = shared_utils.io.fload( fullfile(events_p, un_filename) );
unified_file = shared_utils.io.fload( fullfile(unified_p, un_filename) );

fixation_time = unified_file.opts.STIMULI.setup.go_target.target_duration;

is_go_incorrect = trueat( labels_file.labels, find(labels_file.labels, {'go_trial', 'nogo_choice'}) );

try
  target_event = edf_trials_file.trials('go_target_acquired');
  
  event_times = events_file.event_times;
  event_key = events_file.event_key;
  
  targ_onset_times = event_times(:, event_key('go_target_onset'));
  targ_offset_times = event_times(:, event_key('go_target_offset'));
  targ_acquire_times = event_times(:, event_key('go_target_acquired'));
  
  t0_times = targ_onset_times;
  t1_times = targ_offset_times;
  
%   use_offset_time = is_go_incorrect;
%   use_onset_time = isnan( targ_acquire_times );
%   
%   t0_times = targ_acquire_times;
%   t1_times = targ_acquire_times + (fixation_time * 1e3);
%   
%   t1_times(use_offset_time) = targ_offset_times(use_offset_time);
%   t0_times(use_onset_time) = targ_onset_times(use_onset_time);
  
catch err
  warning( err.message );
  
  continue;
end

t = target_event.time;
t1 = round( t1_times - t0_times );

t1(t1 > max(t)) = max( t );

t0_ind = t == 0;
t1_ind = arrayfun( @(x) find(t == x), t1, 'un', 0 );

if ( nnz(t0_ind) == 0 )
  warning( 'No time 0 found for go target acquired event.' );
  continue;
end

t0_ind = find( t0_ind );

x = target_event.samples('posX');
y = target_event.samples('posY');

go_targ = unified_file.opts.STIMULI.go_target;

targ = go_targ.targets{1};

x_off = targ.x_offset;
y_off = targ.y_offset;
offsets = [ x_off, y_off, x_off, y_off ];
right_bounds = targ.bounds + targ.padding + offsets;

center_x = mean( right_bounds([1, 3]) );
center_y = mean( right_bounds([2, 4]) );

p_ib = rownan( rows(x) );
p_nan = rownan( rows(x) );
c_dist = rownan( rows(x) );
x_dist = rownan( rows(x) );
devs = rownan( rows(x) );

for j = 1:rows(x)
  c_t1_ind = t1_ind{j};
  
  if ( isempty(c_t1_ind) )
    continue; 
  end
  
  select_cols = t0_ind:c_t1_ind;
  
  subset_x = x(j, select_cols);
  subset_y = y(j, select_cols);
  
  is_ib = bfw.bounds.rect( subset_x, subset_y, right_bounds );
  
  p_ib(j) = pnz( is_ib );
  p_nan(j) = pnz( isnan(subset_x) | isnan(subset_y) );
  c_dist(j) = nanmedian( bfw.distance(subset_x, subset_y, center_x, center_y)  );
  x_dist(j) = nanmedian( abs(subset_x - center_x) );
  devs(j) = max( nanstd(subset_x), nanstd(subset_y) );
end

all_inbounds = [ all_inbounds; p_ib ];
all_nan = [ all_nan; p_nan ];
all_dists = [ all_dists; c_dist ];
all_x_dists = [ all_x_dists; x_dist ];
all_devs = [ all_devs; devs ];

append( all_labels, labels_file.labels );

end

assert_ispair( all_inbounds, all_labels );
assert_ispair( all_nan, all_labels );

hwwa.add_data_set_labels( all_labels );

outs.labels = all_labels;
outs.p_inbounds = all_inbounds;
outs.p_nan = all_nan;
outs.dists = all_dists;
outs.xdists = all_x_dists;
outs.devs = all_devs;

end