pro write_xrt_eff_area
; get_xrt_test_effective_area.pro
;
; Updated version to support both single and double filters for test
; file generation - EA.
;
; ==============================================
;
; PROJECT:
;   Solar-B / XRT / XRTpy
;
; NAME:
;   get_xrt_test_effective_area.pro
;
; PURPOSE:
;   Generates effective area test text files to be used for validation
;   comparisons with XRTpy. For each selected XRT channel, it produces a
;   text file containing:
;     - Header metadata (filter, observation date)
;     - Two columns: wavelength [Å] and effective area [cm^2]
;
; CATEGORY:
;   XRTpy testing support
;
; USAGE:
;   1. Modify the `observation_date` variable to your desired test date.
;   2. Run the IDL command:
;        IDL> .compile get_xrt_test_effective_area
;        IDL> write_xrt_eff_area
;
;   Text files will be created in the current directory for each filter
;   listed in the `index` array.
;
; INPUTS:
;   observation_date - (String) Required. Format: 'DD-MMM-YYYY HH:MM:SS'
;                      Example: '22-Sept-2019 22:40:45'
;
; OUTPUTS:
;   Creates text files in the current directory named:
;     <FILTER>_<DATE>_effective_area.txt
;
;   Each file contains:
;     - Filter name
;     - Observation date
;     - Table of wavelength and effective area values
;
; NOTES:
;   - This script uses make_xrt_wave_resp to calculate the wavelength response.
;   - Files are formatted for use by automated test routines in XRTpy.
;   - All filters are processed, including single and supported double-filter combinations.
;   - The filter name in the filename uses dashes instead of slashes (e.g. Al-poly/Ti-poly -> Al-poly-Ti-poly).
;
; CONTACT:
;   Comments, feedback, and bug reports regarding this routine may be
;   directed to this email address: xrt_manager ~at~ head.cfa.harvard.edu
;
; MODIFICATION HISTORY:
;   'v2022-April-25' ;--- (J. Velasquez) Created initial version
;   'v2025-July-14'  ;--- (J. Velasquez) Updated for all filters, added filename safety, clarified usage
;
; ==============================================

  ; === Set the observation date ==============================
  observation_date = '22-Sept-2020 22:40:45'

  ; === Generate wavelength response ===========================
  channel = make_xrt_wave_resp(contam_time=observation_date)

  ; === Create observation date string for filename ============
  year = strmid(observation_date, 0, 2)
  month = strmid(observation_date, 3, 4)
  day = strmid(observation_date, 8, 4)
  observation_date_str = year + month + day ; e.g. 22Sept2019

  ; === List of filter indices (0–13) ===
  index = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]  ; All filters

  ; === Loop through filters and write files ===================
  for i = 0, n_elements(index) - 1 do begin
    filter_name = channel[index[i]].name
    safe_name = str_replace(filter_name, '/', '-')
    filename = safe_name + '_' + observation_date_str + '_effective_area.txt'

    openw, unit, filename, /get_lun

    ; === Write file header ===
    printf, unit, 'Effective Area: IDL Results'
    printf, unit, 'Filter ', filter_name
    printf, unit, 'observation_date ', observation_date
    printf, unit, 'Wavelength (Å)', ' ', 'Effective Area (cm^2)'

    ; === Write wavelength and effective area values ===
    match = where(channel[index[i]].effar.wave ne 0)
    for j = 0, n_elements(match) - 1 do begin
      printf, unit, channel[index[i]].effar.wave[j], ' ', channel[index[i]].effar.eff_area[j]
    endfor

    close, unit, /force
    free_lun, unit
  endfor

end
