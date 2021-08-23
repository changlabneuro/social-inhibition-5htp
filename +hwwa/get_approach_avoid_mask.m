function mask = get_approach_avoid_mask(labels, in_mask, require_initiated)

if ( nargin < 2 || (iscell(in_mask) && isempty(in_mask)) )
  in_mask = rowmask( labels );
end

if ( nargin < 3 )
  require_initiated = true;
end

if ( require_initiated )
  find_initiated = { @find, {'initiated_true'} };
else
  find_initiated = {};
end

mask = fcat.mask( labels, in_mask ...
  , find_initiated{:} ...
  , @findnot, {'tar', '011619'} ...
  , @findnone, 'ephron' ...
  , @findnone, '021819' ... % Eliminate first saline day to match N.
);

mask = union( mask, get_ephron_mask(labels, in_mask) );

end

function mask = get_ephron_mask(labels, mask)

days = {'041019', '041219', '041619', '041819', '042319', '042519' ...
  , '040819', '041119', '041519', '041719', '042219', '042419' ...
  , '042919', '043019', '050219', '050319' ...
};

mask = fcat.mask( labels, mask ...
  , @findor, days ...
  , @find, 'ephron' ...
);

end