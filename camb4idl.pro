;+
;  NAME:
;    camb4idl
;
;  PURPOSE:
;    IDL wrapper for fortran90 CAMB, which you should have installed somewhere
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
done=done[1:n_elements(done)-1]

;MAD Run CAMB, if desired
IF keyword_set(runcamb) THEN BEGIN
   ;MAD Copy over file needed to run CAMB in current directory
   cmd=['cp',filepath('HighLExtrapTemplate_lenspotentialCls.dat', root_dir=getenv('CAMB_DIR')),'.']
   spawn,cmd,/noshell

   ;MAD Run CAMB with new ini file
   cmd=[filepath('camb', root_dir=getenv('CAMB_DIR')),strtrim(outparamfile,2)]
   print,'CAMB4IDL - Calling CAMB with: ' + cmd[0] + ' ' + cmd[1]
   spawn,cmd,camboutput,/noshell

   cambres=parse_cambout(camboutput)
   
   ;MAD Read in outputs
   cambout=read_camb_output(done,outfiles)
   cambres=create_struct(cambout,cambres)
ENDIF

return,cambres
END













