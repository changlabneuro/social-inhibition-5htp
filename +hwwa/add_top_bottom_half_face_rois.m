function [new_rects, new_labs] = add_top_bottom_half_face_rois(rects, labels, mask)

if ( nargin < 3 )
  mask = rowmask( labels );
end

assert_ispair( rects, labels );
face_ind = find( labels, 'face', mask );

new_rects = zeros( numel(face_ind) * 2, 4 );
new_labs = fcat();
stp = 1;

for i = 1:numel(face_ind)
  ind = face_ind(i);
  roi = rects(ind, :);
  h = shared_utils.rect.height( roi );
  
  top_half = [roi(1), roi(2), roi(3), roi(2) + h/2];
  bot_half = [roi(1), roi(2)+h/2, roi(3), roi(4)];
  
  new_rects(stp, :) = top_half;
  append( new_labs, labels, ind );
  setcat( new_labs, 'roi', 'face_top_half', rows(new_labs) );
  
  stp = stp + 1;
  
  new_rects(stp, :) = bot_half;
  append( new_labs, labels, ind );
  setcat( new_labs, 'roi', 'face_bot_half', rows(new_labs) );
  
  stp = stp + 1;
end

end