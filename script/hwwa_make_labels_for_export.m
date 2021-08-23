function results = hwwa_make_labels_for_export(varargin)

defaults = hwwa.get_common_make_defaults();

inputs = { 'meta', 'labels' };
output = 'labels_for_export';

[params, runner] = hwwa.get_params_and_loop_runner( inputs, output, defaults, varargin );

results = runner.run( @main );

end

function out = main(files)

meta_file = shared_utils.general.get( files, 'meta' );
labels_file = shared_utils.general.get( files, 'labels' );

labs = labels_file.labels';

hwwa.add_day_labels( labs );
hwwa.add_data_set_labels( labs );
hwwa.add_drug_labels_by_day( labs );
hwwa.fix_image_category_labels( labs );
hwwa.split_gender_expression( labs );

addcat( labs, 'monkey' );
setcat( labs, 'monkey', lower(meta_file.monkey) );

out = struct();
out.unified_filename = meta_file.unified_filename;
[out.labels, out.categories] = categorical( labs );

end