function write_anova_tables(measure, grps, p, prefix)

grp_vars = cellfun( @(x) measure(:, x), grps, 'un', false );

[~, tbl, stats] = anovan( measure.data, grp_vars ...
  , 'model', 'full' ...
  , 'varnames', grps ...
  );

prefix2 = strjoin( incat(measure, grps), '_' );

tbl = table( tbl );

[c, m, ~, g] = multcompare( stats, 'dimension', 1:numel(grps) );

names = arrayfun( @(x) g(x), c(:, 1:2) );
c_cell = arrayfun( @(x) {x}, c );
c_cell(:, 1:2) = names;
sig = c(:, end) < 0.05;
sig_names = c_cell(sig, :);

c_tbl = table( c_cell );

fname1 = fullfile( p, sprintf('%s_%s_anova.csv', prefix2, prefix) );
fname2 = fullfile( p, sprintf('%s_%s_posthoc.csv', prefix2, prefix) );

writetable( tbl, fname1 );
writetable( c_tbl, fname2 );

end