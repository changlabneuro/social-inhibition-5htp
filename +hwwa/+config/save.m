
function save(conf)

%   SAVE -- Save the config file.

hwwa.util.assertions.assert__is_config( conf );
const = hwwa.config.constants();
fprintf( '\n Config file saved\n\n' );
save( fullfile(const.config_folder, const.config_filename), 'conf' );

hwwa.config.load( '-clear' );

end