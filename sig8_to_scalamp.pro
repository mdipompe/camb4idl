;+
;  NAME:
;    sig8_to_scalamp
;
;  PURPOSE:
;    CAMB doesn't take sigma8 as an imput, instead uses scalar
;    amplitude.  This will tell you what scalar amp to use to give you
;    desired sigma8 at z=0
;
;  USE:
;    scalamp=sig8_to_scalamp(sigma8,paramfile='default_params.ini')
;
;  INPUT:
;    sigma8 - desired sigma8 at z=0
;
;  OPTIONAL INPUT:
;    paramfile - CAMB parameter file (defaults to default_params.ini)
;    h0 - little h (H0/100), defaults to 0.702
;    omega_m - Omega_matter, defaults to 0.273
;    omega_l - Omega_lambda, defaults to 0.727
;    omega_b - omega_baryon, defaults to 0.046
;    spec_ind - power spectrum spectral index, defaults to 0.96
;
;  OUTPUT:
;    scalamp - the value you should set for scalar amplitude in CAMB
;              to get desired sigma_8
;
;  NOTES:
;    Relies on CAMB and CAMB4IDL wrapper.
;
;  HISTORY:
;    11-7-15 - Written - MAD (Dartmouth)
;    11-30-15 - Added power spectrum spectral index n - MAD (Dartmouth)
;-
FUNCTION sig8_to_scalamp,sigma8,paramfile=paramfile,$
                           h0=h0,omega_m=omega_m,omega_l=omega_l,$
                           omega_b=omega_b,spec_ind=spec_ind

;MAD Set defaults  
IF ~keyword_set(paramfile) THEN paramfile='default_params.ini'
IF ~keyword_set(h0) THEN h0=0.702
IF ~keyword_set(omega_m) THEN omega_m=0.275
IF ~keyword_set(omega_b) THEN omega_b=0.046
IF ~keyword_set(omega_l) THEN omega_l=0.725
IF ~keyword_set(spec_ind) THEN spec_ind=0.96
omega_dm=omega_m-omega_b

;MAD CALL CAMB with test scalar_amp 2e-9
res=camb4idl(/runcamb, paramfile=paramfile,output_root='thesearecambtestfiles_', $
             outparamfile='thesearecambtestfiles_out.ini', $
             get_scalar='T',get_transfer='T', $
             get_tensor='F',do_lensing='T', $
             do_nonlinear=0, $
             scalar_amp=2e-9, $
             scalar_spectral_index=spec_ind, $
             transfer_redshift=0, $
             transfer_filename='transfer_0.dat',$
             transfer_matterpower='transfer_0.dat', $
             use_physical='F', $
             hubble=h0*100., $
             omega_baryon=omega_b, $
             omega_cdm=omega_dm, $
             omega_lambda=omega_l)

;MAD Parse output and get measured sigma_8
currentsig=-99
i=0
WHILE (currentsig EQ -99) DO BEGIN
   xx=strsplit(res.raw_output[i],' ',/extract)
   IF (n_elements(xx) GT 6) THEN BEGIN
      IF (xx[0] EQ 'at' AND xx[1] EQ 'z' AND xx[2] EQ '=' $
          AND xx[3] EQ '0.000' AND xx[4] EQ 'sigma8' $
          AND xx[5] EQ '(all') THEN currentsig=float(xx[8])
   ENDIF
   IF (i EQ 99) THEN currentsig=-1
   i=i+1
ENDWHILE

;MAD Exit if couldn't parse the output properly
IF (currentsig EQ -1) THEN BEGIN
   print,'SIG8_TO_SCALAMP: Couldn''t find sigma_8 in test CAMB'
   print,'                 output...exiting, returning -1'
   return,'-1'
ENDIF

;MAD Rescale scalar amp to get desired sigma_8
scalaramp=(2e-9)*((sigma8/currentsig)^2)

;MAD Delete the files test CAMB run created
files=file_search('thesearecambtestfiles_*')
cmd=['rm',files,'HighLExtrapTemplate_lenspotentialCls.dat']
spawn,cmd,/noshell

return,scalaramp
END
