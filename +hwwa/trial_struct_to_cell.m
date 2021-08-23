function [cells, map] = trial_struct_to_cell(data, map)

%   TRIAL_STRUCT_TO_CELL -- Convert struct trial data to a cell matrix.

if ( nargin == 1 )
  map = containers.Map();
end

fields = fieldnames( data );
rows = numel( data );
n_fields = numel( fields );

error_prefix = 'error__';
is_error_field = cellfun( @(x) ~isempty(strfind(x, error_prefix)), fields );
n_error_fields = sum( is_error_field );

cols = n_fields - n_error_fields;

cells = cell( rows, cols );

assigned_cols = [];

for i = 1:rows
  
  trial = data(i);
  next_col = 1;
  
  for j = 1:n_fields
  
    field = fields{j};

    assign_field = field;
    assign_value = trial.(field);

    %   convert 'error__broke_fixation' to 'broke_fixation'
    if ( is_error_field(j) )
      assign_field = 'errors';
      
      if ( assign_value )
        assign_value = field(numel(error_prefix)+1:end);
      else
        assign_value = 'no_errors';
      end
    end

    if ( isKey(map, assign_field) )
      assign_col = map(assign_field);
    else
      assert( ~any(assigned_cols == next_col), 'Already assigned to column %d.', next_col );
      assign_col = next_col;
      map(assign_field) = next_col;
      assigned_cols(end+1) = next_col;
      next_col = next_col + 1;
    end

    cells( i, assign_col ) = { assign_value };
  end
end

end