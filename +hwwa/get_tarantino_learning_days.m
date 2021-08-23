function days = get_tarantino_learning_days()

days = { '0117', '0118', '0121', '0122', '0123' };
days = cellfun( @(x) [x, '19'], days, 'un', 0 );

end