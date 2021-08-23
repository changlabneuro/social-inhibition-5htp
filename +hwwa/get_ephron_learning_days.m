function days = get_ephron_learning_days()

days = { '0318', '0319', '0320', '0321', '0322' };
days = cellfun( @(x) [x, '19'], days, 'un', 0 );

end