;+
;  NAME:
;    read_camb_output
;
;  PURPOSE:
;    Read in files that CAMB writes out
;
;  USE:
;    result=read_camb_output(used,file_struct)
;
;  INPUT:
;    used - string from what_did_camb_do listing what was calculated
;    file_struct - structure containing possible file names
;    
;  OUTPUT:
;    Returns a structure with all the CAMB data
;
;  HISTORY:
;    5-22-15 - Written - MAD (UWyo)
;-
FUNCTION read_camb_output,used,files
  tmp=create_struct('tmp',0)
  ;Reorder what was used into scaler, vector, tensor, transfer, lensing
  xx=where(used EQ 'T')
  

  IF (used[0] EQ 'S') THEN BEGIN     ;Scalar cls
     scalcls=read_matrix(files.scalar_output_file)
     tmp=create_struct(tmp,'scalar_cls',scalcls)
  ENDIF
  IF (used[1] EQ 'V') THEN BEGIN     ;Vector cls
     vectcls=read_matrix(files.vector_output_file)
     tmp=create_struct(tmp,'vector_cls',vectcls)
  ENDIF
  IF (used[2] EQ 'T') THEN BEGIN     ;Tensor cls
     tenscls=read_matrix(files.tensor_output_file)
     tmp=create_struct(tmp,'tensor_cls',tenscls)
  ENDIF
  IF (used[0] EQ 'S' AND used[2] EQ 'T') THEN BEGIN     ;combined cls
     totcls=read_matrix(files.total_output_file)
     tmp=create_struct(tmp,'total_cls',totcls)
  ENDIF
  IF (used[4] EQ 'l') THEN BEGIN     ;lensed cls
     lenscls=read_matrix(files.lensed_output_file)
     tmp=create_struct(tmp,'lensed_cls',lenscls)
  ENDIF
  IF (used[0] EQ 'S' AND used[2] EQ 'T' $      ;combined lensed cls
      AND used[4] EQ 'l') THEN BEGIN 
     lenstotcls=read_matrix(files.lensed_total_output_file)
     tmp=create_struct(tmp,'lensed_total_cls',lenstotcls)
  ENDIF
  IF (used[3] EQ 'tr') THEN BEGIN    ;transfer functions
     FOR i=0L,n_elements(files.transfer_filename)-1 DO BEGIN
        tf=read_matrix(files.transfer_filename[i])
        tagname='transfer_func_' + strtrim(i+1,2)
        tmp=create_struct(tmp,tagname,tf)
     ENDFOR
     FOR i=0L,n_elements(files.transfer_matterpower)-1 DO BEGIN
        ps=read_matrix(files.transfer_matterpower[i])
        tagname='matter_power_' + strtrim(i+1,2)
        tmp=create_struct(tmp,tagname,ps)
     ENDFOR
  ENDIF

  remove_tags,tmp,'tmp',out

  return,out
END
