function labs = get_unit_labels(unit)

labs = fcat.with( {'channel', 'id', 'ms_channel', 'region'} );
setcat( labs, 'channel', unit.channel );
setcat( labs, 'id', unit.id );

if ( isfield(unit, 'ms_channel') )
  setcat( labs, 'ms_channel', unit.ms_channel );
end

if ( isfield(unit, 'region') )
  setcat( labs, 'region', unit.region );
end

end