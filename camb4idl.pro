;These are all the functions that camb4idl depends on, camb4idl itself
;is at the bottom...
;----------------------------------------------------------------------------------------
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



;+
;  NAME:
;    read_matrix
;  PURPOSE:
;    read in a matrix from a text file
;
;  USE:
;    matrix = read_square_matrix('filename.txt')           
;
;  INPUT:
;    file - string name of file to read
;
;  RETURNS:
;    matrix - the array read from the file
;
;  HISTORY:
;    5-12-15 - Written - MAD (UWyo)
;-
FUNCTION read_matrix,file

  ;Get dimensions by counting number of lines
  openr,1,file
  WHILE (not EOF(1)) DO BEGIN
     line=''
     readf,1,line
     xx=strsplit(line,' ',/extract)
     IF (n_elements(C) EQ 0) THEN C=double(xx) ELSE C=[[C],[double(xx)]]
  ENDWHILE
  close,1

  return,C
END




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





;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME:
;    REMOVE_TAGS
;       
; PURPOSE:
;    remove the specified tags from input structure
;
; CALLING SEQUENCE:
;    remove_tags, oldstruct, tagnames, newstruct
;
; INPUTS: 
;    oldstruct: the original structure
;    tagnames: the names of tags to be removed (can be an array)
;
; OPTIONAL INPUTS:
;    NONE.
;
; KEYWORD PARAMETERS:
;    NONE.
;       
; OUTPUTS: 
;    newstruct: the new structure without tags.
;
; OPTIONAL OUTPUTS:
;    NONE
;
; CALLED ROUTINES:
;    
; 
; PROCEDURE: 
;    
;	
;
; REVISION HISTORY:
;    ????? Judith Racusin
;    25-OCT-2000 Modified to handle arbitrary tag types. Also error 
;          handling. Erin Scott Sheldon
;       
;                                      
;-                                       
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO remove_tags, struct, tagnames, newstruct

  IF n_params() EQ 0 THEN BEGIN 
      print,'Syntax - remove_tags, oldstruct, tagnames, newstruct'
      print
      print,'Use doc_library,"remove_tags"  for more help.'  
      return
  END

  ;; Figure out which tags get removed

  tags=tag_names(struct)
  n=n_elements(tags)
  tagnames=strupcase(tagnames)
  nt=n_elements(tagnames)
  IF nt EQ 1 THEN BEGIN
      t=where(tags NE tagnames[0],nw) 
      IF nw EQ n THEN BEGIN
          print,'-----------------------------------------------------'
          message,'Tag did not match, structure unchanged',/inf
          print,'-----------------------------------------------------'
          newstruct = struct
          return
      ENDIF 
  ENDIF ELSE BEGIN 
      match,tags,tagnames,m
      IF m[0] EQ -1 THEN BEGIN
          print,'-------------------------------------------------'
          message,'No tags matched, structure unchanged',/inf
          print,'-------------------------------------------------'
          newstruct=struct
          return
      ENDIF 
      nm=n_elements(m)
      IF nm EQ n THEN BEGIN 
          print,'-------------------------------------------------------------'
          message,'This would remove all tags! structure unchanged',/inf
          print,'-------------------------------------------------------------'
          newstruct=struct
          return
      ENDIF 
      t=lindgen(n)
      remove, m, t
  ENDELSE 
      
  ;; create new structure
  tags=tags[t]
  n=n_elements(tags)

  newstruct=create_struct(tags[0],struct[0].(t[0]))
  
  FOR i=1L, n-1 DO newstruct = create_struct(temporary(newstruct), $
                                             tags[i], struct[0].(t[i]) )

  newstruct=replicate( temporary(newstruct), n_elements(struct) )
  struct_assign,struct,newstruct

  return
END







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
;----------------------------------------------------------------------------------------




;+
;  NAME:
;    camb4idl
;
;  PURPOSE:
;    IDL wrapper for fortran90 CAMB, which you should have installed
;    somewhere (the env variable CAMB_DIR should point there)
;
;  USE:
;    result=camb4idl(paramfile='input_params.ini',outparamfile='output_params.ini',/runcamb,[camb_keys...])
;
;  INPUT:
;    
;  OPTIONAL INPUT:
;    paramfile - an input CAMB parameter file.  Defaults to
;                'default_params.ini'
;    outparamfile - the paramfile that will be written out and passed
;                   to camb.  Identical to input if no keywords set.
;    output_root, etc. - inputs to CAMB.  See paramfile.  A small
;                        minority not currently accepted (see
;                        commented lines in function def below)
;
;  KEYWORDS:
;    runcamb - Will call CAMB and run on output param file.  Must have
;              environment variable CAMB_DIR defined to location of CAMB executable.
;
;  OUTPUT:
;    IF runcamb is set, outputs all of the CAMB files as set in param
;    file. Returns a structure with all of this information, and more.
;
;  NOTES: 
;    Param file should have spaces before/after '=' signs (this is
;    required for accurate parsing)
;
;    Does not currently support multiple scalar amplitudes, spectral
;    indices, etc.
;
;    DOES support an array of redshifts for the transfer
;    function/matter power spectrum, and the code will handle naming
;    them and writing them to files (in fact this functionality is
;    largely why I wrote it!).  
;
;    You shouldn't put a list of redshifts/files in the default param
;    file by hand (even if just two) - let the code do it for you.
;
;    Note that if you set
;    transfer_redshift to an array with n_elements > 1, you must also
;    set transfer_filename and transfer_matterpower to string arrays
;    of the same length.  CAMB will only do up to 150 redshifts at a time.
;
;  HISTORY:
;    5-22-15 - Written - MAD (UWyo)
;-
FUNCTION camb4idl,paramfile=paramfile,outparamfile=outparamfile,runcamb=runcamb,$
                 output_root=output_root, $
                 get_scalar_cls=get_scalar_cls, $
                 get_vector_cls=get_vector_cls, $
                 get_tensor_cls=get_tensor_cls, $
                 get_transfer=get_transfer, $
                 do_lensing=do_lensing, $
                 do_nonlinear=do_nonlinear, $
                 l_max_scalar=l_max_scalar, $
                 l_max_tensor=l_max_tensor, $
;                 k_eta_max_tensor=k_eta_max_tensor, $
                 use_physical=use_physical, $
                 ombh2=ombh2, $
                 omch2=omch2, $
                 omnuh2=omnuh2, $
                 omk=omk, $
                 hubble=hubble, $
                 w=w, $
                 cs2_lam=cs2_lam, $
                 omega_baryon=omega_baryon, $
                 omega_cdm=omega_cdm, $
                 omega_lambda=omega_lambda, $
                 omega_neutrino=omega_neutrino, $
                 temp_cmb=temp_cmb, $
                 helium_fracion=helium_fracion, $
                 massless_neutrinos=massless_neutrinos, $
                 nu_mass_eigenstates=nu_mass_eigenstates, $
                 massive_neutrinos=massive_neutrinos, $
                 share_delta_neff=nare_delta_neff, $
                 nu_mass_fractions=nu_mass_fractions, $
                 ss_neutrinos=ss_neutrions, $
                 nu_mass_degeneracies=nu_mass_degeneracies, $
                 initial_power_num=initial_power_num, $
                 pivot_scalar=pivot_scalar, $
                 pivot_tensor=pivot_tensor, $
                 scalar_amp=scalar_amp, $
                 scalar_spectral_index=scalar_spectral_index, $
                 scalar_nrun=scalar_nrun, $
                 scalar_nrunrun=scalar_nrunrun, $
                 tensor_spectral_index=tensor_spectral_index, $
                 tensor_nrun=tensor_nrun, $
                 tensor_parameterization=tensor_parameterization, $
                 initial_ratio=initial_ratio, $
                 tensor_amp=tensor_amp, $
                 reionization=reionization, $
                 re_use_optical_depth=re_use_optical_depth, $
                 re_optical_depth=re_optical_depth, $
                 re_redshift=re_redshift, $
                 re_delta_redshift=re_delta_redshift, $
                 re_ionization_frac=re_ionization_frac, $
                 RECFAST_fudge=RECFAST_fudge, $
                 RECFAST_fudge_He=RECFAST_fudge_He, $
                 RECFAST_Heswitch=RECFAST_heswitch, $
                 RECFAST_Hswitch=RECFAST_Hswitch, $
                 cosmorec_runmode=cosmorec_runmode, $
                 cosmorec_accuracy=cosmorec_accuracy, $
                 cosmorec_fdm=cosmorec_fdm, $
                 initial_condition=initial_condition, $
                 initial_vector=initial_vector, $
                 vector_mode=vector_mode, $
                 COBE_normalize=COBE_normalize, $
                 CMB_outputscale=CMB_outputscale, $
                 transfer_high_precision=transfer_high_precision, $
                 transfer_kmax=transfer_kmax, $
                 transfer_k_per_logint=transfer_k_per_logint, $
                 transfer_num_redshifts=transfer_num_redshifts, $
                 transfer_interp_matterpower=transfer_interp_matterpower, $
                 transfer_redshift=transfer_redshift, $
                 transfer_filename=transfer_filename, $
                 transfer_matterpower=transfer_matterpower, $
                 transfer_power_var=transfer_power_var, $
                 scalar_output_file=scalar_output_file, $
                 vector_output_file=vector_output_file, $
                 tensor_output_file=tensor_output_file, $
                 total_output_file=total_output_file, $
                 lensed_output_file=lensed_output_file, $
                 lensed_total_output_file=lensed_total_output_file, $
                 lens_potential_output_file=lens_potential_output_file, $
                 FITS_filename=FITS_filename, $
;                 do_lensing_bispectrum=do_lensing_bispectrum, $
                 do_primordial_bispectrum=do_primordial_bispectrum, $
                 bispectrum_nfields=bispectrum_nfields, $
                 bispectrum_slice_base_L=bispectrum_slice_base_L, $
;                 bispectrum_ndelta=bispectrum_ndelta, $
;                 bispectrum_delta=bispectrum_delta, $
                 bispectrum_do_fisher=bispectrum_do_fisher, $
                 bispectrum_fisher_noise=bispectrum_fisher_noise, $
                 bispectrum_fisher_noise_pol=bispectrum_fisher_noise_pol, $
                 bispectrum_fisher_fwhm_arcmin=bispectrum_fisher_fwhm_arcmin, $
                 bispectrum_full_output_file=bispectrum_full_output_file, $
                 bispectrum_full_output_sparse=bispectrum_full_output_sparse, $
                 bispectrum_export_alpha_beta=bispectrum_export_alpha_beta, $
                 feedback_level=feedback_level, $
                 derived_parameters=derived_parameters, $
                 lensing_method=lensing_method, $
                 accurate_BB=accurate_BB, $
                 massive_nu_approx=massive_nu_approx, $
                 accurate_polarization=accurate_polarization, $
                 accurate_reionization=accurate_reionization, $
                 do_tensor_neutrinos=do_tensor_neutrinos, $
                 do_late_rad_truncation=do_late_rad_truncation, $
                 halofit_version=halofit_version, $
                 number_of_threads=number_of_threads, $
                 high_accuracy_default=high_accuracy_default, $
                 accuracy_boost=accuracy_boost, $
                 l_accuracy_boost=l_accuracy_boost, $
                 l_sample_boost=l_sample_boost

;MAD Set default input and output param files
IF ~keyword_set(paramfile) THEN paramfile='default_params.ini'
IF ~keyword_set(outparamfile) THEN outparamfile='my_params.ini'

;MAD Read in paramfile as defaults
defaults=read_defaults(paramfile)

;MAD Get list of keywords to check
params=routine_info('camb4idl',/parameters,/functions)
params=params.kw_args

;MAD Initialize output param file
openw,1,outparamfile

;MAD Done is a bitnumber indicating what CAMB did
done=' '

;Loop over lines, keeping defaults or resetting as needed
FOR i=0L,n_elements(defaults)-1 DO BEGIN
   ;MAD check if line includes keyword that should be reset
   result=get_keyword(defaults[i],params)
   IF (result EQ 'notkey') THEN BEGIN   ;Line is a comment or other, just print as-is
      printf,1,defaults[i]
   ENDIF ELSE BEGIN
      ;MAD Get output root from file or call
      IF (result EQ 'output_root') THEN BEGIN
         IF ~keyword_set(output_root) THEN outroot=get_value(defaults[i]) ELSE $
            outroot=output_root
      ENDIF

      ;MAD Get number of redshifts from file or call
      IF (result EQ 'transfer_num_redshifts') THEN BEGIN
         IF ~keyword_set(transfer_num_redshifts) THEN nz=double(get_value(defaults[i])) ELSE $
            nz=transfer_num_redshifts
      ENDIF

      ;MAD Some keywords require special action, check if this is one
      special=check_special(result)

      ;MAD Set a flag if the line is a
      ;potential keyword and has been set in call
      cmd=execute('IF (n_elements(' + strtrim(result,2) + ') NE 0) THEN flag=1')
      IF (n_elements(flag) EQ 0) THEN BEGIN   ;Not in call, should stay as default
         printf,1,defaults[i]                

         IF (special EQ 3 OR special EQ 5) THEN BEGIN
            filestruc=get_outfile(result,get_value(defaults[i]),outroot)
            IF (n_elements(outfiles) EQ 0) THEN outfiles=filestruc ELSE $
               outfiles=create_struct(outfiles,filestruc)
         ENDIF

         ;MAD Add to 'done' what the defaults say to do
         IF (special EQ 4) THEN BEGIN
            did=what_did_camb_do(result,get_value(defaults[i]))
            done=[done,did]
         ENDIF

      ENDIF ELSE BEGIN       ;Set in call, new value should be used
         ;MAD if value should be changed, set tmp to new value for writing
         x=execute('tmp='+strtrim(result,2))

         IF (special EQ 3 OR special EQ 5) THEN BEGIN
            filestruc=get_outfile(result,tmp,outroot)
            IF (n_elements(outfiles) EQ 0) THEN outfiles=filestruc ELSE $
               outfiles=create_struct(outfiles,filestruc)
         ENDIF


         ;MAD Add to 'done' what the call says to do
         IF (special EQ 4) THEN BEGIN
            did=what_did_camb_do(result,tmp)
            done=[done,did]
         ENDIF

         ;MAD Some keywords need a number following them, add if needed
         IF (special NE 1 AND special NE 2 AND special NE 5) THEN BEGIN ;Doesn't need modification, write
            printf,1,strtrim(result,2) + '  =  ' + strtrim(tmp,2)     
         ENDIF ELSE IF (special EQ 1) THEN BEGIN         ;Needs number added to keyword in param file
            printf,1,strtrim(result,2) + '(1)  =  ' + strtrim(tmp,2)  
         ENDIF ELSE IF (special EQ 2 OR special EQ 5) THEN BEGIN         ;Is a redshift line, loop over them
            FOR j=0L,n_elements(transfer_redshift)-1 DO BEGIN    
               string=strtrim(result,2) + '(' + strtrim(j+1,2) + $
                      ')    = ' + strtrim(tmp[j],2)
               printf,1,string
            ENDFOR
         ENDIF
      ENDELSE
      undefine,flag
   ENDELSE
ENDFOR
close,1

;MAD Reorder what was done into scalar, vector, tensor, trans, lensing 
done=done[1:n_elements(done)-1]
done2=['S','V','T','tr','l']
FOR i=0L,n_elements(done2)-1 DO BEGIN
   xx=where(done2[i] EQ done)
   IF (xx[0] EQ -1) THEN done2[i]='N'
ENDFOR
done=done2

;MAD If nothing is set to be calculated, then you must have just
;wanted to make a parameter file, and set runcamb by mistake
IF (done[0] EQ 'N' AND done[1] EQ 'N' AND done[2] EQ 'N') THEN oops=1 ELSE oops=0
IF (keyword_set(runcamb) AND oops EQ 1) THEN BEGIN
   print, 'CAMB4IDL - You said to run CAMB but don''t want to calculate anything? Returning...'
   return,oops
ENDIF

;MAD Run CAMB, if desired
IF (keyword_set(runcamb) AND oops EQ 0) THEN BEGIN
   ;MAD Copy over file needed to run CAMB in current directory
   cmd=['cp',filepath('HighLExtrapTemplate_lenspotentialCls.dat', root_dir=getenv('CAMB_DIR')),'.']
   spawn,cmd,/noshell

   ;MAD Run CAMB with new ini file
   cmd=[filepath('camb', root_dir=getenv('CAMB_DIR')),strtrim(outparamfile,2)]
   print,'CAMB4IDL - Calling CAMB with: ' + cmd[0] + ' ' + cmd[1]
   spawn,cmd,camboutput,/noshell

   cambres=parse_cambout(camboutput,done,nz)
   
   ;MAD Read in outputs
   cambout=read_camb_output(done,outfiles)
   cambres=create_struct(cambout,cambres)

   return,cambres
ENDIF

END





