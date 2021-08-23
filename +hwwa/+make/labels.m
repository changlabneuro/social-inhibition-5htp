function labels_file = labels(files)

%   LABELS -- Create labels file.
%
%     file = ... labels( files ) creates the labels intermediate file from 
%     the unified intermediate file contained in `files`. `files` is a 
%     containers.Map or struct with an entry for 'unified'.
%
%     The output file contains an fcat object that represents trial data in
%     a human readable way; conceptually, it is simply a categorical array,
%     and can be converted to one with the categorical() function.
%
%     See also fcat
%
%     IN:
%       - `files` (containers.Map, struct)
%     FILES:
%       - 'unified'
%     OUT:
%       - `events_file` (struct)

assert( ~isempty(which('fcat')), ['This function depends on the categorical' ...
  , ' repository available at: https://github.com/nfagan/categorical'] );

hwwa.validatefiles( files, 'unified' );

unified_file = shared_utils.general.get( files, 'unified' );
unified_filename = hwwa.try_get_unified_filename( unified_file );

reprocess_fields = { 'cue_delay', 'trial_number', 'trial_outcome', 'trial_type' };
ignore_fields = [ {'reaction_time', 'events', 'target_displacement', 'reward'}, reprocess_fields ];
ignore_fields{end+1} = 'reward_size_cue_file';

data = unified_file.DATA;

str_delays = num_field_to_str( data, 'cue_delay', 'delay__' );
str_trials = num_field_to_str( data, 'trial_number', 'trial__' );

all_fields = fieldnames( data );
use_fields = setdiff( all_fields, ignore_fields );

f = addcat( fcat(), [use_fields(:)', reprocess_fields] );

for j = 1:numel(use_fields)
  field = use_fields{j};
  
  values = { data.(field) };
  
  if ( strcmp(field, 'target_type') )
    values(strcmp(values, 'social')) = { 'social_target' };
  end
  
  setcat( f, field, values );
end

setcat( f, 'cue_delay', str_delays );
setcat( f, 'trial_number', str_trials );

trial_outs = { data(:).trial_outcome };
trial_types = { data(:).trial_type };

empty_outs = cellfun( @isempty, trial_outs );
trial_outs(empty_outs) = { 'no' };

trial_outs = cellfun( @(x) [x, '_choice'], trial_outs, 'un', false );
trial_types = cellfun( @(x) [x, '_trial'], trial_types, 'un', false );

setcat( f, 'trial_outcome', trial_outs );
setcat( f, 'trial_type', trial_types );

%
%   date info
%

addcat( f, {'date', 'drug', 'unified_filename'} );
setcat( f, 'date', unified_file.date );
setcat( f, 'unified_filename', unified_filename );

%
%   derived labels
%

hwwa.add_correct_labels( f );
hwwa.add_prev_trial_correct_labels( f );
hwwa.add_initiated_labels( f );
hwwa.add_drug_labels( f );
hwwa.add_trial_type_change_labels( f );

labels_file = struct();
labels_file.labels = f;
labels_file.unified_filename = unified_filename;

end

function out = num_field_to_str(data, fieldname, prefix)
out = arrayfun( @(x) [prefix, num2str(x.(fieldname))], data, 'un', false );
end