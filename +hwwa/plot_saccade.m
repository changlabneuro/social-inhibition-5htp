function axs = plot_saccade(t, x, y, starts, stops)

ax1 = subplot( 1, 2, 1 );
cla( ax1 );
plot( ax1, t, x );
title( ax1, 'x' );

for i = 1:numel(starts)
  hold( ax1, 'on' );
  shared_utils.plot.add_vertical_lines( ax1, t(starts(i)) );
  shared_utils.plot.add_vertical_lines( ax1, t(stops(i)), 'r--' );
end

ax2 = subplot( 1, 2, 2 );
cla( ax2 );
plot( ax2, t, y );
title( ax2, 'y' );

for i = 1:numel(starts)
  hold( ax2, 'on' );
  shared_utils.plot.add_vertical_lines( ax2, t(starts(i)) );
  shared_utils.plot.add_vertical_lines( ax2, t(stops(i)), 'r--' );
end

axs = [ ax1, ax2 ];

end