function varargout = social_divide_scrambled(data, labels, each, mask)

if ( nargin < 4 )
  mask = rowmask( labels );
end

a = 'not-scrambled';
b = 'scrambled';

opfunc = @rdivide;
sfunc = @(x) nanmean( x, 1 );

[varargout{1:nargout}] = dsp3.sbop( data, labels, each, a, b, opfunc, sfunc, mask );

if ( nargout > 1 )
  setcat( varargout{2}, 'scrambled_type', sprintf('%s - %s', a, b) );
end

end