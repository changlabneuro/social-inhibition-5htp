function days = get_hitch_learning_days()

days = { '0116', '0117', '0118', '0121', '0122' };
days = cellfun( @(x) [x, '19'], days, 'un', 0 );

end