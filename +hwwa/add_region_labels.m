function labs = add_region_labels(labs, prefix)

un_p = hwwa.get_intermediate_dir( 'unified' );

C = combs( labs, 'unified_filename' );

reg_cat = 'region';

addcat( labs, reg_cat );

for i = 1:numel(C)
  
  un_f = C{i};
  
  un_fullfile = fullfile( un_p, un_f );
  
  if ( ~shared_utils.io.fexists(un_fullfile) )
    fprintf( '\n Unified file "%s" not found.', un_fullfile );
    continue;
  end
  
  un_file = shared_utils.io.fload( un_fullfile );
  
  if ( ~isfield(un_file, 'region_map') || isempty(un_file.region_map) )
    fprintf( '\n Unified file "%s" does not have a region map.', un_f );
    continue;
  end
  
  reg_map = un_file.region_map;
  
  regs = fieldnames( reg_map );
  
  for j = 1:numel(regs)
    region = regs{j};
    channels = arrayfun( @(x) channel_to_str(prefix, x), reg_map.(region), 'un', false );
    
    for k = 1:numel(channels)
      ind = find( labs, channels{k} );
      if ( isempty(ind) )
        fprintf( '\n No channels matched "%s" for "%s".', channels{k}, un_f );
        continue; 
      end
      setcat( labs, reg_cat, region, ind );
    end
  end
end

end

function str = channel_to_str(prefix, n)

if ( n < 10 )
  str = sprintf( '%s0%d', prefix, n );
else
  str = sprintf( '%s%d', prefix, n );
end

end