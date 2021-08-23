function labs = add_day_labels( labs )

[I, C] = findall( labs, 'date' );

addcat( labs, 'day' );

for i = 1:numel(I)
  
  try
    s = datestr( C{i}, 'mmddyy' );
  catch err
    s = 'day_NaN';
  end
  
  setcat( labs, 'day', s, I{i} );
end

end