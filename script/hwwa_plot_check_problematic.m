conf = hwwa.config.load();
conf.PATHS.data_root = '/Volumes/My Passport/NICK/Chang Lab 2016/hww_gng/data';

use_files = { 'test_', '1022', '1023', '1024', '1025', '1003', '1001', '1011', '1012' };
% use_files = { 'test_' };

% use_files =  { '1025' };

hwwa.make_edf_trials( 'files_containing', use_files ...
  , 'config',       conf ...
  , 'event',        {'go_target_onset'} ...
  , 'look_back',    0 ...
  , 'look_ahead',   1e3 ...
  , 'overwrite',    true ...
  , 'append',       true ...
);

%%
%   , 'files_containing',   'test_human_' ...
outs = hwwa_check_problematic_go_trials_looped( ...
    'config',             hwwa_set_local_data_root() ...
  , 'kind',               'stats' ...
  , 'files_containing',   'test' ...
  , 'is_parallel',        true ...
  , 'target_event',       'go_target_onset' ...
);

%%

% data_col = outs.data_key('met_fixation_criterion');
data_col = outs.data_key('blink_duration');

pltdat = outs.data(:, data_col);
pltlabs = outs.labels';
traces = outs.traces;

un_filenames = combs( pltlabs, 'unified_filename' );
test_filenames = shared_utils.cell.containing( un_filenames, 'test' );

addsetcat( pltlabs, 'data-set', 'monk' );
setcat( pltlabs, 'data-set', 'human', find(pltlabs, test_filenames) );

%%

pl = plotlabeled.make_common();

mask = fcat.mask( pltlabs ...
  , @findnone, 'no_choice' ...
  , @find, {'correct_false', 'go_trial'} ...
);

xcats = { 'trial_type' };
% gcats = { 'data-set' };
gcats = { 'unified_filename' };
pcats = { 'trial_outcome', 'correct' };

pl.bar( pltdat(mask), pltlabs(mask), xcats, gcats, pcats );

%%

pltdat = outs.data(:, outs.data_key('met_fixation_criterion'));

mask = fcat.mask( pltlabs ...
  , @findnone, 'no_choice' ...
  , @find, {'correct_false', 'go_trial'} ...
);

[y, I] = keepeach( pltlabs', {'data-set'}, mask );
p = rowzeros( numel(I) );

for i = 1:numel(I)
  
  subset = pltdat(I{i});
  
  nans = isnan( subset );
  is1 = subset == 1;
  
  N = nnz( ~nans );
  p(i) = sum( is1 ) / N;
  
end

bar( p );