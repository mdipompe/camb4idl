;+
;  NAME:
;    get_keyword
;
;  PURPOSE:
;    Splits up a string from the parameter file to pull out a
;    potential keyword.  Identifies as a keyword or not.
;
;  USE:
;    result=get_keyword(string,list)
;
;  INPUT:
;    string - the string to parse
;    list - the list of possible keywords to check against
;    
;  OUTPUT:
;    Returns string - 'notkey' or the name of the keyword.
;
;  HISTORY:
;    5-22-15 - Written - MAD (UWyo)
;-
FUNCTION get_keyword,line,list
  split1=strsplit(line,' ',/extract)
  split2=strsplit(split1[0],'(',/extract)
  match=where(list EQ STRUPCASE(split2[0]))
  IF (match[0] EQ -1) THEN out='notkey' ELSE out=split2[0]
  return,out
END
