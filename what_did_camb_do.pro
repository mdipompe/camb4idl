;+
;  NAME:
;    what_did_camb_do
;
;  PURPOSE:
;    Determines what your CAMB run did, based on 'T' and 'F'
;    statements for things like get_scalar_cls, get_tensor_cls,
;    do_lensing, etc.
;
;  USE:
;    result=what_did_camb_do(key,val)
;
;  INPUT:
;    key - the keyword you're checking
;    val - the value of the keyword (T or F)
;    
;  OUTPUT:
;    Returns 5 element string (scalar, vector, tensor, transfer,
;    lensing), set to N if not done, and S/V/T/tr/l if done.
;
;  HISTORY:
;    5-22-15 - Written - MAD (UWyo)
;-
FUNCTION what_did_camb_do,key,val
  to_do=['get_scalar_cls', $
         'get_vector_cls', $
         'get_tensor_cls', $
         'get_transfer', $
         'do_lensing']

  xx=where(key EQ to_do)
  IF (xx EQ 0) THEN out='S'
  IF (xx EQ 1) THEN out='V'
  IF (xx EQ 2) THEN out='T'
  IF (xx EQ 3) THEN out='tr'
  IF (xx EQ 4) THEN out='l'

  IF (val EQ 'F') THEN out='N'

  return,out
END
