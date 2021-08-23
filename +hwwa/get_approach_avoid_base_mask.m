function mask = get_approach_avoid_base_mask(labels, mask_func)

require_initiated = true;

mask = hwwa.get_approach_avoid_mask( labels, {}, require_initiated );
mask = mask_func( labels, mask );

end