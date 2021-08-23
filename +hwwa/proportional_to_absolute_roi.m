function rect = proportional_to_absolute_roi(image_roi, proportional_roi)

im_width = shared_utils.rect.width( image_roi );
im_height = shared_utils.rect.height( image_roi );

min_x = im_width * proportional_roi(1) + image_roi(1);
min_y = im_height * proportional_roi(2) + image_roi(2);

abs_width = im_width * shared_utils.rect.width( proportional_roi );
abs_height = im_height * shared_utils.rect.height( proportional_roi );

rect = [ min_x, min_y, min_x + abs_width, min_y + abs_height ];

end