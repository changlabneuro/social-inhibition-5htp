function roi_file = rois(files)

unified_file = shared_utils.general.get( files, 'unified' );

roi_file = struct();
roi_file.unified_filename = unified_file.unified_filename;

go_targ = go_target_roi( unified_file );
nogo_cue = nogo_cue_roi( unified_file );

roi_file.nogo_cue = nogo_cue;
roi_file.nogo_cue_padded = pad_roi( nogo_cue, nogo_padding(unified_file) );
roi_file.go_target = go_targ;
roi_file.go_target_padded = pad_roi( go_targ, go_target_padding(unified_file) );
roi_file.go_target_displacement = go_target_displacement( unified_file );

end

function rect = go_target_displacement(unified_file)

displacement = unified_file.opts.STIMULI.setup.go_target.displacement;

if ( isempty(displacement) )
  displacement = 0;
end

if ( numel(displacement) == 1 )
  displacement(end+1) = displacement(1);
end

rect = [ 0, 0, displacement(1), displacement(2) ];

end

function pad = go_target_padding(unified_file)

pad = unified_file.opts.STIMULI.setup.go_target.target_padding;

end

function pad = nogo_padding(unified_file)

pad = unified_file.opts.STIMULI.setup.nogo_cue.target_padding;

end

function roi = pad_roi(roi, xy)

x = xy(1);

if ( numel(xy) == 1 )
  y = x;
else
  y = xy(2);
end

cx = mean( roi([1, 3]) );
cy = mean( roi([2, 4]) );

w = (roi(3) - roi(1)) + x*2;  % padding is +/- x
h = (roi(4) - roi(2)) + y*2;  % padding is +/- y

roi = [ cx - w/2, cy - h/2, cx + w/2, cy + h/2 ];

end

function roi = go_target_roi(unified_file)

roi = unified_file.opts.STIMULI.go_target.vertices;

end

function roi = nogo_cue_roi(unified_file)

roi = unified_file.opts.STIMULI.nogo_cue.vertices;

end