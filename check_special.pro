;+
;  NAME:
;    check_special
;
;  PURPOSE:
;    Identifies keywords that require special action in the main code.    
;
;  USE:
;    result=check_special(keyword)
;
;  INPUT:
;    keyword - the (string) name of the keyword to check against the list
;    
;  OUTPUT:
;    Returns 0 (not special), 1 (needs additional string), 2 (related to redshift
;    list), 3 (contains name of output file), 4 (contains info on
;    what calculations CAMB did), or 5 (both redshift and output file)
;
;  HISTORY:
;    5-22-15 - Written - MAD (UWyo)
;-
FUNCTION check_special,keyword
  mod_list=['scalar_amp', $
            'scalar_spectral_index', $
            'scalar_nrun', $
            'scalar_nrunrun', $
            'tensor_spectral_index', $
            'tensor_nrun', $
            'initial_ratio', $
            'tensor_amp']

  redshift_list=['transfer_redshift']
 

  output_list=['scalar_output_file', $
               'vector_output_file', $
               'tensor_output_file', $
               'total_output_file', $
               'lensed_output_file', $
               'lensed_total_output_file', $
               'lens_potential_output_file']

  do_list=['get_scalar_cls', $
           'get_vector_cls', $
           'get_tensor_cls', $
           'get_transfer', $
           'do_lensing']

  z_and_out_list=['transfer_filename', $
                 'transfer_matterpower']

  out=0
  xx=where(keyword EQ mod_list)
  IF (xx[0] NE -1) THEN out = 1 
  xx=where(keyword EQ redshift_list)
  IF (xx[0] NE -1) THEN out = 2
  xx=where(keyword EQ output_list)
  IF (xx[0] NE -1) THEN out = 3
  xx=where(keyword EQ do_list)
  IF (xx[0] NE -1) THEN out = 4
  xx=where(keyword EQ z_and_out_list)
  IF (xx[0] NE -1) THEN out = 5
  
  return,out
END
