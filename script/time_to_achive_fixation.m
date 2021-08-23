un_p = hwwa.get_intermediate_dir( 'unified' );
lab_p = hwwa.get_intermediate_dir( 'labels' );
evt_p = hwwa.get_intermediate_dir( 'events' );
lab_mats = hwwa.require_intermediate_mats( lab_p );
lab_files = cellfun( @shared_utils.io.fload, lab_mats );
labs = arrayfun( @(x) x.labels, lab_files, 'un', false );
labs = extend( labs{1}', labs{2:end} );

evt_files = arrayfun( @(x) fullfile(evt_p, x.unified_filename), lab_files, 'un', false );
evts = cellfun( @shared_utils.io.fload, evt_files );
events = cell2mat( {evts(:).event_times}' );
evt_key = evts(1).event_key;

%%

time_diffs = [];

for i = 1:numel(evt_files)
  evts = shared_utils.io.fload( evt_files{i} );
  labs = shared_utils.io.fload( fullfile(lab_p, evts.unified_filename) );
  unified = shared_utils.io.fload( fullfile(un_p, evts.unified_filename) );
  
  if ( i == 1 )
    all_labs = fcat.like( labs.labels );
  end
  
  event_times = evts.event_times;
  cue_times = event_times(:, cue_ind);
  fix_times = event_times(:, fix_ind);
  N = numel( fix_times );
  
  fix_dur = unified.opts.TIMINGS.fixations.fix_square;
  
  fix_ind = evts.event_key('fixation_on');
  cue_ind = evts.event_key('go_nogo_cue_onset');
  
  stp = 1;
  next = 1;
  
  time_diff = [];
  
  while ( next <= N )
    
    c_fix = fix_times(stp);
    c_cue = cue_times(stp);
    
    while ( isnan(c_cue) && next < N )
      next = next + 1;
      c_cue = cue_times(next);
    end
    
    time_diff = [ time_diff; (c_cue - c_fix - fix_dur) ];
    
    stp = next + 1;
    next = stp;
  end  
  
  labs = one( labs.labels' );
  repeat( labs, numel(time_diff)-1 );
  
  time_diffs = [ time_diffs; time_diff ];
  append( all_labs, labs );
  
%   time_diff = evts.event_times(:, cue_ind) - evts.event_times(:, fix_ind);
%   time_diff = time_diff - fix_dur;
%   
%   time_diffs = [ time_diffs; time_diff ];
%   append( all_labs, labs.labels );
end

%%

dat = labeled( time_diffs, all_labs );

% only( dat, 'initiated_true' );

means = eachindex( dat', {'initiated', 'date'}, @rownanmean );
devs = each( dat', {'initiated', 'date'}, @(x) std(x, [], 1) );

pl = plotlabeled();
pl.error_func = @plotlabeled.nansem;

pl.bar( means, 'drug', 'trial_outcome', 'trial_outcome' );


