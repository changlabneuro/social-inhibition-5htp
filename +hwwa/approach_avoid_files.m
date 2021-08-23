function files = approach_avoid_files()

files = [ hwwa.get_image_5htp_days(), hwwa.get_image_saline_days() ];
files = cellstr( hwwa.to_date(files) );

end