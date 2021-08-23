function m = make_mask(l, f)
m = f( l, rowmask(l) );
end