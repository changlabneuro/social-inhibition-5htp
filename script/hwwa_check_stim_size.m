stim_sessions = bfw.get_sessions_by_stim_type( [], 'cache', true );

%%

labels_p = hwwa.gid( 'labels' );
unified_p = hwwa.gid( 'unified' );

un_mats = shared_utils.io.findmat( unified_p );
un_mats = shared_utils.io.filter_files( un_mats, {'11-F', '13-F', '14-F'} );

image_files_per_day = cell( size(un_mats) );

sz_per_day = zeros( numel(un_mats), 2 );
verts_per_day = zeros( numel(un_mats), 4 );
labs = fcat();
did_persist = false( numel(un_mats), 1 );

xs = zeros( numel(un_mats), 1 );
ys = zeros( size(xs) );

for i = 1:numel(un_mats)
  unified_file = shared_utils.io.fload( un_mats{i} );
  labels_file = shared_utils.io.fload( fullfile(labels_p, unified_file.unified_filename) );
  
  verts = unified_file.opts.STIMULI.go_target.vertices;
  
  sz_per_day(i, :) = [verts(3)-verts(1), verts(4)-verts(2)];
  verts_per_day(i, :) = verts;
  xs(i) = mean( [verts(3), verts(1)] );
  ys(i) = mean( [verts(4), verts(2)] );
  did_persist(i) = unified_file.opts.STRUCTURE.persist_images_through_reward;
  
  addcat( labels_file.labels, 'monkey' );
  setcat( labels_file.labels, 'monkey', lower(unified_file.opts.META.monkey) );
  
  append1( labs, labels_file.labels );
  
%   image_files_per_day{i} = unique( {unified_file.DATA.image_file} );
end

%%