function [h, r, p] = add_corr(ax, x, y, varargin)

defaults = struct();
defaults.alpha = 0.05;
defaults.add_text = true;
defaults.text_prefix = '';
defaults.p_format = '%0.3f';
defaults.r_format = '%0.2f';

params = shared_utils.general.parsestruct( defaults, varargin );

r_format = use_default_if_empty( params, 'r_format', defaults );
p_format = use_default_if_empty( params, 'p_format', defaults );

[r, p] = corr( x, y, 'rows', 'complete' );

xlims = get( ax, 'xlim' );
ylims = get( ax, 'ylim' );
xticks = get( ax, 'xtick' );

ps = polyfit( x, y, 1 );
y = polyval( ps, xticks );

cstate = get( ax, 'nextplot' );
set( ax, 'nextplot', 'add' );
h = plot( ax, xticks, y );
set( ax, 'nextplot', cstate );

h.Annotation.LegendInformation.IconDisplayStyle = 'off';

coord_func = @(x) ((x(2)-x(1)) * 0.75) + x(1);

xc = coord_func( xlims );
yc = y(end);

if ( yc < ylims(1) )
  yc = ylims(1);
elseif ( yc > ylims(2) )
  yc = ylims(2);
end

p_txt = sprintf( p_format, p );
r_txt = sprintf( r_format, r );
txt = sprintf( '%sR = %s, p = %s', params.text_prefix, r_txt, p_txt );

if ( p < params.alpha ), txt = sprintf( '%s *', txt ); end

if ( params.add_text )
  t = text( ax, xc, yc, txt );
  set( t, 'color', get(h, 'color') );
end

end

function v = use_default_if_empty(params, field, defaults)

v = params.(field);
if ( isempty(v) )
  v = defaults.(field);
end

end