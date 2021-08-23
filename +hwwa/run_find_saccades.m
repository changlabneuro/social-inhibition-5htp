function start_stops = run_find_saccades(x_deg, y_deg)

smooth_func = @(data) smoothdata( data, 'smoothingfactor', 0.05 );

vel_thresh = 100;
dur_thresh = 50;

start_stops = hwwa.find_saccades( x_deg, y_deg, 1e3, vel_thresh, dur_thresh, smooth_func );

end