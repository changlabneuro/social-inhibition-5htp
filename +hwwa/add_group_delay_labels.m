function labs = add_group_delay_labels( labs, starts, stops )

assert( numel(starts) == numel(stops), 'Number of starts must match number of stops.' );

delay_strs = combs( labs, 'cue_delay' );
delays = shared_utils.container.cat_parse_double( 'delay__', delay_strs );

group_cat = 'gcue_delay';
group_labs = cell( size(labs, 1), 1 );

leftovers = true( size(group_labs) );

for i = 1:numel(starts)
  
  start = starts(i);
  stop = stops(i);
  
  group_ind = delays >= start & delays <= stop;
  matching_ind = find( labs, delay_strs(group_ind) );
  
  group_lab = sprintf( 'grouped_delay__%s-%s', num2str(start), num2str(stop) );
      
  group_labs(matching_ind) = { group_lab };
  leftovers(matching_ind) = false;
end

group_labs(leftovers) = { 'grouped_delay__NaN' };

addcat( labs, group_cat );
setcat( labs, group_cat, group_labs );

end