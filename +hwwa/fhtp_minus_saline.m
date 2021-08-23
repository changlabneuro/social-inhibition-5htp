function varargout = fhtp_minus_saline(data, labels, each, mask)

if ( nargin < 4 )
  mask = rowmask( labels );
end

a = '5-htp';
b = 'saline';

opfunc = @minus;
sfunc = @(x) nanmean( x, 1 );

[varargout{1:nargout}] = dsp3.sbop( data, labels, each, a, b, opfunc, sfunc, mask );

if ( nargout > 1 )
  setcat( varargout{2}, 'drug', sprintf('%s - %s', a, b) );
end

end