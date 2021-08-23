function write_anova_summary_tables( measure, grps, save_p )

means = eachindex( measure', grps, @rowmean );
devs = each( measure', grps, @(x) std(x, [], 1) );

means_tbl = table( Container.from(means), grps );
devs_tbl = table( Container.from(devs), grps );

fname = strjoin( incat(means, grps), '_' );
fname = strrep( fname, '.', '_' );

writetable( means_tbl, fullfile(save_p, sprintf('means_%s.csv', fname)), 'WriteRowNames', true );
writetable( devs_tbl, fullfile(save_p, sprintf('devs_%s.csv', fname)), 'WriteRowNames', true );

end