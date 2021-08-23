function meta_file = meta(files)

import shared_utils.struct.field_or;

hwwa.validatefiles( files, 'unified' );

unified_file = shared_utils.general.get( files, 'unified' );

meta_file = struct();
meta_file.unified_filename = unified_file.unified_filename;
meta_file.monkey = field_or( unified_file.opts.META, 'monkey', '' );

end