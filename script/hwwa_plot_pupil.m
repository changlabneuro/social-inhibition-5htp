function hwwa_plot_pupil(pupil, labels, varargin)

assert_ispair( pupil, labels );

defaults = hwwa.get_common_plot_defaults( hwwa.get_common_make_defaults() );
defaults.mask_func = @(labels, mask) mask;

params = hwwa.parsestruct( defaults, varargin );

mask = get_base_mask( labels, params.mask_func );

[norm_pupil, norm_labels] = saline_normalize( pupil, labels', mask );

compare_to_one( norm_pupil, norm_labels, 'pupil_size', params );
change_from_saline( norm_pupil, norm_labels, params );

end

function [normed, norm_labels] = saline_normalize(pupil, labels, mask)

norm_each = { 'monkey' };
[normed, norm_labels] = hwwa.saline_normalize( pupil, labels, norm_each, mask );

end

function compare_to_one(norm_pupil, norm_labels, subdir, params)

sr_outs = dsp3.signrank1( norm_pupil, norm_labels', {} ...
  , 'signrank_inputs', {1} ...  % test against median of 1.
);

if ( params.do_save )
  data_p = get_data_path( params, 'analyses', subdir );
  dsp3.save_signrank1_outputs( sr_outs, data_p, 'drug' );
end

end

function change_from_saline(normed, norm_labels, params)

fcats = {};
gcats = { 'drug' };
pcats = {};

pl = plotlabeled.make_common();
axs = pl.boxplot( normed, norm_labels, gcats, pcats );

subdir = 'pupil_size';

if ( params.do_save )
  save_p = get_data_path( params, 'plots', subdir );
  dsp3.req_savefig( gcf, save_p, norm_labels, [gcats, pcats, fcats] );
end

end

function mask = get_base_mask(labels, mask_func)

mask = hwwa.get_approach_avoid_mask( labels );
mask = mask_func( labels, mask );

end

function p = get_data_path(params, data_type, varargin)

p = hwwa.approach_avoid_data_path( params, data_type, varargin{:} );

end

function p = get_plot_path(params, varargin)

p = get_data_path( params, 'plots', varargin{:} );

end