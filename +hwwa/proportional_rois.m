function rois = proportional_rois(rois)

% From keynote
image_width = 768;
image_height = 768;

for i = 1:size(rois, 1)
  width = shared_utils.rect.width( rois(i, :) );
  height = shared_utils.rect.height( rois(i, :) );
  
  off_x = rois(i, 1);
  off_y = rois(i, 2);
  
  frac_width = width / image_width;
  frac_height = height / image_height;
  frac_x = off_x / image_width;
  frac_y = off_y / image_height;
  
  rois(i, :) = [ frac_x, frac_y, frac_x + frac_width, frac_y + frac_height ];
end

end