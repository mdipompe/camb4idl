;+
;  NAME:
;    parse_cambout
;
;  PURPOSE:
;    Take (some) stuff CAMB prints to terminal and parse into a
;    structure. Most just goes in as a string array to be parsed by
;    the user.
;
;  USE:
;    result=parse_cambout(input)
;
;  INPUT:
;    input - the stuff CAMB prints to screen
;    done - string from what_did_camb_do indicating what was
;           calculated and thus returned
;    numz - the number of zs for the transfer function/power spectrum    
;
;  OUTPUT:
;    Returns structure with tags for some CAMB output and values
;
;  HISTORY:
;    5-22-15 - Written - MAD (UWyo)
;-
FUNCTION parse_cambout,input,done,numz
  taglist=['z_reion', $
           'Age_universe', $
           'z_star']
  tagval=dblarr(n_elements(taglist))

  j=0
  FOR i=0L,n_elements(input)-1 DO BEGIN
     IF (i EQ 0 OR i EQ 11 OR i EQ 12) THEN BEGIN
        tmp=strsplit(input[i],' ',/extract)
        tagval[j]=double(tmp[n_elements(tmp)-1])
        j=j+1
     ENDIF
  ENDFOR
  
  FOR i=0L,n_elements(taglist)-1 DO BEGIN
     IF (i EQ 0) THEN out=create_struct(taglist[i],tagval[i]) ELSE $
        out=create_struct(out,taglist[i],tagval[i])
  ENDFOR

  IF (done[3] EQ 'tr') THEN BEGIN
     IF (done[0] EQ 'S' AND done[2] EQ 'T') THEN i=26
     IF (done[0] EQ 'S' AND done[2] EQ 'N') THEN i=25
     IF (done[0] EQ 'N' AND done[2] EQ 'T') THEN i=26
     IF (done[0] EQ 'N' AND done[2] EQ 'N') THEN i=25

     sigma_8=dblarr(numz)
     j=0
     WHILE (j LT n_elements(sigma_8)) DO BEGIN
        tmp=strsplit(input[i],' ',/extract)
        sigma_8[j]=double(tmp[n_elements(tmp)-1])
        i=i+1
        j=j+1
     ENDWHILE
     out=create_struct(out,'sigma_8',sigma_8)
  ENDIF

  out=create_struct(out,'raw_output',input)

  return,out
END
