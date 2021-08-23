function out = unified_filename_to_datetime(filename)

if ( ischar(filename) )
  out = convert( filename );
else
  out = cellfun( @convert, filename, 'un', 0 );
end

out = datetime( out );

end

function converted = convert(s)

converted = s(1:end-4);

dash_ind = find( converted == '-' );
if ( numel(dash_ind) ~= 4 )
  error( 'Invalid format for string: "%s".', s );
end

converted(dash_ind(3:4)) = ':';

end