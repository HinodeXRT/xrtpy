pro write_xrt_tr
; ==============================================
;
; PROJECT:
;   Solar-B / XRT / XRTpy
;
; NAME:
;
;   MAKE_XRT_TR
;
; PURPOSE:
;   Generates temperature response test text files to be used for validation
;   comparisons with XRTpy. For each selected XRT channel, it produces a
;   text file containing:
;     - Header metadata (filter, observation date, abundance)
;     - Two columns: temperature [K] and temperature response [DN cm^5 pix^-1 s^-1]
;
; CATEGORY:
;   XRTpy testing support
;
; CALLING SEQUENCE:
;   get_xrt_test_temperature_response
;
;
; INPUTS:
;   None directly. You must set `observation_date` and optional `abundance_model`
;   manually inside the script for now.
;
; MODIFIABLE PARAMETERS:
;   observation_date - (String) Required. Format: 'DD-MMM-YYYY HH:MM:SS'
;                      Example: '22-Sept-2023 21:45:45'
;
;   abundance_model  - (String) One of 'photospheric', 'hybrid', or 'coronal'
;
;   chn_filename     - GENX filename to use (e.g. 'xrt_channels_v0017')
;
; KEYWORDS:
;
;   Index - This keyword is a list of numbers that correspond to
;           a filter. A text file will be created for the
;           filters in the list. Currently set to all filters.
;           Al-mesh = 0
;           Al-poly = 1
;           C-poly = 2
;           Ti-poly = 3
;           Be-thin = 4
;           Be-med = 5
;           Al-med = 6
;           Al-thick = 7
;           Be-thick = 8
;
;
; OUTPUTS:
;   Creates text files in the current directory named:
;     <FILTER>_<DATE>_temperature_response_<ABUNDANCE>.txt
;
; NOTES:
;   - Make sure CHIANTI abundance flags match the chosen model.
;   - This script uses make_xrt_wave_resp and make_xrt_temp_resp under the hood.
;   - Files are formatted for use by automated test routines in XRTpy.
;
;
; INFORMATION EXTENSION:
;
;   make_xrt_wave_resp - Reference make_xrt_temp_resp.pro for the
;                        procedure to calculate the effective areas and
;                        spectral responses for a set of XRT x-ray channel
;                        accounting for some thickness of the CCD
;                        contamination layer. The spectral response is
;                        directly calculated from the effective area, and both are
;                        functions of wavelength.
;
;   make_xrt_temp_resp - Reference make_xrt_temp_resp.pro for the
;                        procedure to calculate the temperature response.
;
; CONTACT:
;
;   Comments, feedback, and bug reports regarding this routine may be
;   directed to this email address:
;   xrt_manager ~at~ head.cfa.harvard.edu 
;   For comments, feedback, and bug reports regarding xrtpy realted items -
;   directed to this email address: xrtpy ~at~ cfa.harvard.edu
;
; MODIFICATION HISTORY:
;
;   'v2025-April-1' ;--- (J.Velasquez) Create
;                     get_xrt_test_temperature_response to create
;                     testing text files of the temperature response
;                     for XRTpy
;   NOTE: After first run:
;         IDL> retall
;         IDL> .compile get_xrt_test_temperature_response
;         IDL> write_xrt_tr
;
; ==============================================
;
;
; === Insert an observation date to test ======================

  observation_date = '22-Sept-2023 21:45:45'

  ; ==== Individual string of the observation date created for text file title==
  year = strmid(observation_date, 0, 2)
  month = strmid(observation_date, 3, 4)
  day = strmid(observation_date, 8, 4)
  observation_date_str = year + month + day ; format for string file title

  ;Generate wavelength response using selected channel file and contamination date
  wave_resp = make_xrt_wave_resp(contam_time=observation_date, chn_filename='xrt_channels_v0017')

  ;Calculate temperature response with selected CHIANTI model
  ;'photospheric': temp_resp = make_xrt_temp_resp(wave_resp, /chianti_defaul, /photospheric)
  ;'hybrid':       temp_resp = make_xrt_temp_resp(wave_resp, /chianti_defaul, /hybrid)
  ;'coronal':      temp_resp = make_xrt_temp_resp(wave_resp, /chianti_defaul, /coronal)
  temp_resp = make_xrt_temp_resp(wave_resp,/chianti_defaul,/photospheric)
  
  ; === index is a list of numbers corresponding to each filter ===
  index = [0, 1, 2, 3, 4, 5, 6, 7, 8]

  ; === IDL for loop creating a text file containing what is stated in
  ; OUTPUTS ===
  for i = 0, n_elements(index) - 1 do begin
    openw, unit, wave_resp[index[i]].name + '_' + observation_date_str + '_temperature_response_photospheric.txt', /get_lun

    printf, unit, 'Temperature Response: IDL Results'
    printf, unit , 'Abundance_model: Photospheric'
    printf, unit, 'Filter ', wave_resp[index[i]].name
    printf, unit, 'observation_date ', observation_date
    printf, unit, 'Temperature', ' ', 'Temperature Response'

    for j = 0, n_elements(temp_resp[index[i]].temp_resp) - 1 do begin ; -1
      printf, unit, temp_resp[index[i]].temp[j], ' ', temp_resp[index[i]].temp_resp[j]
    endfor

    close, unit, /force
    free_lun, unit
  endfor

end
