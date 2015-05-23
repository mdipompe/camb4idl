;+
;  NAME:
;    get_value
;
;  PURPOSE:
;    Some cases need to save the value of a parameter, even if not set
;    in call (i.e. need the default).  This gets it.
;
;  USE:
;    result=get_value(line)
;
;  INPUT:
;    line - the parameter file line to read
;    
;  OUTPUT:
;    Returns value of parameter as a string
;
;  HISTORY:
;    5-22-15 - Written - MAD (UWyo)
;-
FUNCTION get_value,line
  split=strsplit(line,' ',/extract)
  out=split[n_elements(split)-1]
  return,out
END
