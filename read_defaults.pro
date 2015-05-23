;+
;  NAME:
;    read_defaults
;
;  PURPOSE:
;    Read in initial parameter file, split into string array
;
;  USE:
;    result=read_defaults(paramfile)
;
;  INPUT:
;    file - string name of parameter file
;    
;  OUTPUT:
;    Returns string array, each element a line of the parameter file
;
;  HISTORY:
;    5-22-15 - Written - MAD (UWyo)
;-
FUNCTION read_defaults,file
  openr,1,file
  line=''
  WHILE ~EOF(1) DO BEGIN
     readf,1,line
     IF (n_elements(lines) EQ 0) THEN lines=line ELSE $
        lines=[lines,line]
  ENDWHILE
  close,1  

  return,lines
END
