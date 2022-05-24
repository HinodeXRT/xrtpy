pro write_xrt_eff_area
; ==============================================
;
; PROJECT:
;       Solar-B / XRT / XRTpy 
;
; NAME:
;       
;       MAKE_XRT_EFF_AREA
;
; CATEGORY:
;       
;       XRTpy 
;
; PURPOSE:
;
;       Produces testing formated text files containing header
;       information i.e Filter,observation_data, and two rows
;       of data; wavelength and effective area for a set of XRT
;       x-ray channels paired with thicknesses of the CCD 
;       contamination layer.
;
; INPUTS:
;
;       INPUT1   - [Mandatory] observation_date, (string)
;		      String formated observation date is required
;		      to calculate the effective area. Format of
;		      observation date 'DD-MMMM-YYYY HH:MM:SS'.
;
; KEYWORDS:
;
;        Index -  This keyword is a list of number that coressponse to
;	       	        a filter. A text file will be created for the
;			filters in the list. Currently set to all filters.
;	 		Al-mesh = 0
;	 		Al-poly = 1
;	 		C-poly = 2
;	 		Ti-poly = 3
;	 		Be-thin = 4
;	 		Be-med = 5
;	 		Al-med = 6
;	 		Al-thick = 7
;	 		Be-thick = 8
;
; OUTPUTS:
;
;       Return - The text files contain header information i.e Filter,
;              	      observation_data, and two rows of data; wavelength and
; 		      effective area. The format of creating these text files is
;    		      important for a proper testing run.
;                     -		Effective Area (cm**2)
;		      -		Wavelength (Angstroms)
; INFORMATION EXTENSION:
;
;      make_xrt_wave_resp  - Reference make_xrt_temp_resp.pro for the
;			   procedure to calcuate the effective areas and spectral responses
;			   for a set of XRT x-ray channel accounting for some thickness
;       		   of the CCD contamination layer. The spectral response is
;       		   directly calculated from the effective area, and both are
;       		   functions of wavelength.
;
; CONTACT:
;
;       Comments, feedback, and bug reports regarding this routine may be
;       directed to this email address:
;                xrt_manager ~at~ head.cfa.harvard.edu
;
; MODIFICATION HISTORY:
;
;       'v2022-April-25' ;--- (J.Velasquez) Create
;                             get_xrt_test_effective_area to create
;                             testing text files of the effective area
;                             for XRTpy
;
;  ==============================================
;
;
; === Insert an observation date to test ======================

  observation_date ='22-Sept-2019 22:40:45'


  channel = make_xrt_wave_resp(contam_time = observation_date)


; ==== Individual string of the observation date create for text file
; title ==

  year = strmid( observation_date, 0,2)
  month = strmid( observation_date, 3,4)
  day = strmid( observation_date, 8,4)
  observation_date_str = year+month+day

; === index is a list numbers corresponding to each filter ===
  index = [0,1,2,3,4,5,6,7,8]

; === IDL for loop creating a text file containing what is stated in
; OUTPUTS ===
  
  for i=0,n_elements(index) -1 do begin
     openw,unit,channel[index[i]].name+'_'+ observation_date_str +'_effective_area.txt',/get_lun

     printf, unit, 'Filter ',channel[index[i]].name
     printf, unit, 'observation_date ', observation_date_str 
     printf, unit, 'wavelength', ' ', 'effective area '
     match = where(channel[index[i]].effar.wave ne 0)
     
     for j=0,n_elements(match)-1 do begin
        printf, unit, channel[index[i]].effar.wave[j],' ',channel[index[i]].effar.eff_area[j]
     endfor
     
     close, unit, /force

  endfor

end
