pro write_xrt_tr
; ===============================================================================
; get_xrt_test_temperature_response.pro
;
; PROJECT:
;   Solar-B / Hinode / XRT / XRTpy Testing Utilities
;
; NAME:
;   write_xrt_tr
;
; PURPOSE:
;   Generates temperature response test text files for XRTpy validation.
;   Supports both single and double filter channels. For each filter, creates a text file with:
;     - Header metadata (filter, observation date, abundance model)
;     - Two columns: temperature [K] and temperature response [DN cm^5 pix^-1 s^-1]
;
; CALLING SEQUENCE:
;   IDL> .compile get_xrt_test_temperature_response
;   IDL> write_xrt_tr
;
; REQUIRED SETUP:
;   - `make_xrt_wave_resp` and `make_xrt_temp_resp` must be in your IDL path
;   - Ensure CHIANTI environment is properly configured
;
; PARAMETERS:
;   observation_date - String, e.g., '22-Sept-2010 21:45:45'
;   abundance        - One of 'photospheric', 'hybrid', or 'coronal'
;
; OUTPUT:
;   Creates a series of files named:
;     <FILTER>_<DATE>_temperature_response_<ABUNDANCE>.txt
;   in the current working directory.
;
; NOTE:
;   - Double filters use a hyphen in filenames (e.g., Al-poly/Ti-poly → Al-poly-Ti-poly)
;   - Text files are formatted for automated comparison testing against XRTpy outputs
;
; CONTACT:
;   xrtpy ~at~ cfa.harvard.edu
;
; MODIFICATION HISTORY:
;   v2025-July-14 — J.Velasquez — Updated to support abundance selection and better documentation
; ==============================================================================

  ; === Insert an observation date to test ======================
  observation_date = '22-Sept-2010 21:45:45'

  ; === Create observation date string for text file naming ====
  year = strmid(observation_date, 0, 2)
  month = strmid(observation_date, 3, 4)
  day = strmid(observation_date, 8, 4)
  observation_date_str = year + month + day ; e.g. 22Sept2010

  ; === Generate wavelength response ===
  wave_resp = make_xrt_wave_resp(contam_time=observation_date, chn_filename='xrt_channels_v0017')

  ; === Define abundance model as variable ===
  abundance = 'hybrid'

  ; === Choose CHIANTI model ===
  if abundance eq 'photospheric' then $
      temp_resp = make_xrt_temp_resp(wave_resp, /chianti_default, /photospheric) $
  else if abundance eq 'hybrid' then $
      temp_resp = make_xrt_temp_resp(wave_resp, /chianti_default, /hybrid) $
  else if abundance eq 'coronal' then $
      temp_resp = make_xrt_temp_resp(wave_resp, /chianti_default, /coronal) $
  else begin
      message, 'Invalid abundance model. Use "photospheric", "hybrid", or "coronal".'
  endelse

  ; === List of filter indices (0–13) ===
  index = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]  ; All filters

  ; === Loop through filters and write files ===
  for i = 0, n_elements(index) - 1 do begin
    filter_name = wave_resp[index[i]].name
    safe_name = str_replace(filter_name, '/', '-')
    filename = safe_name + '_' + observation_date_str + '_temperature_response_' + abundance + '.txt'

    openw, unit, filename, /get_lun

    printf, unit, 'Temperature Response: IDL Results'
    printf, unit, 'Abundance_model: ' + abundance
    printf, unit, 'Filter ', filter_name
    printf, unit, 'observation_date ', observation_date
    printf, unit, 'Temperature', ' ', 'Temperature Response'

    for j = 0, n_elements(temp_resp[index[i]].temp_resp) - 1 do begin
      printf, unit, temp_resp[index[i]].temp[j], ' ', temp_resp[index[i]].temp_resp[j]
    endfor

    close, unit, /force
    free_lun, unit
  endfor

end
