;+
;  NAME:
;    get_outfile
;
;  PURPOSE:
;    builds a structure of output file names for use in read in parser
;
;  USE:
;    result=get_outfile(key,val,root)
;
;  INPUT:
;    key - the keyword of the output file
;    val - the name of the output file
;    root - the root of the output files
;    
;  OUTPUT:
;    Returns structure with all the possible out file names (not all
;    necesarily used, depends on your CAMB settings)
;
;  HISTORY:
;    5-22-15 - Written - MAD (UWyo)
;-
FUNCTION get_outfile,key,val,root
  output_list=['transfer_filename', $
               'transfer_matterpower', $
               'scalar_output_file', $
               'vector_output_file', $
               'tensor_output_file', $
               'total_output_file', $
               'lensed_output_file', $
               'lensed_total_output_file', $
               'lens_potential_output_file']

  xx=where(key EQ output_list)
  
  return,create_struct(output_list[xx],root+'_'+val)
END
