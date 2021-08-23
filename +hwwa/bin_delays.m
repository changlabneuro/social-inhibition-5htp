function [delay_starts, delay_stops] = bin_delays(delays, n_delays)

n_split = floor( numel(delays) / n_delays );

delay_starts = zeros( n_delays, 1 );
delay_stops = zeros( size(delay_starts) );

stp = 1;

for i = 1:n_delays
  delay_starts(i) = delays(stp);
  
  if ( i < n_delays )
    M = min( stp+n_split-1, numel(delays) );
  else
    M = numel(delays);
  end
  
  delay_stops(i) = delays(M);
  
  stp = stp + n_split;
end

end