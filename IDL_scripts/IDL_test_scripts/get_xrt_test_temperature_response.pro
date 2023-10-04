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
; CATEGORY:
;
;   XRTpy
;
; PURPOSE:
;
;   Produces temperature response testing formatted text files
;   containing header information i.e Filter, observation_data, and two rows
;   of data; temperature and temperature response for a set of XRT
;   x-ray channels.
;
; INPUTS:
;
;   INPUT1 - [Mandatory] observation_date, (string)
;            String formatted observation date is required
;            to calculate the temperature response. Format of
;            observation date: 'DD-MMMM-YYYY HH:MM:SS'.
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
; OUTPUTS:
;
;   Return - The text files contain header information i.e Filter,
;            observation_data, and two rows of data; temperature and
;            temperature response. The format of creating these text files is
;            important for a proper testing run.
;            - Temperature (degrees Kelvin, K)
;            - Temperature Response (el cm^5 s^-1 pix^-1)
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
;
; MODIFICATION HISTORY:
;
;   'v2022-April-25' ;--- (J.Velasquez) Create
;                     get_xrt_test_temperature_response to create
;                     testing text files of the temperature response
;                     for XRTpy
;
; ==============================================
;
;
; === Insert an observation date to test ======================

  observation_date = '24-Oct-2023 21:25:45'

  ; ==== Individual string of the observation date created for text file title==
  year = strmid(observation_date, 0, 2)
  month = strmid(observation_date, 3, 4)
  day = strmid(observation_date, 8, 4)
  observation_date_str = year + month + day ; format for string file title

  wave_resp = make_xrt_wave_resp(contam_time=observation_date)
  temp_resp = make_xrt_temp_resp(wave_resp, /chianti)

  ; === index is a list of numbers corresponding to each filter ===
  index = [0, 1, 2, 3, 4, 5, 6, 7, 8]

  ; === IDL for loop creating a text file containing what is stated in
  ; OUTPUTS ===
  for i = 0, n_elements(index) - 1 do begin
    openw, unit, wave_resp[index[i]].name + '_' + observation_date_str + '_temperature_response.txt', /get_lun

    printf, unit, 'Temperature Response: IDL Results'
    printf, unit, 'Filter ', wave_resp[index[i]].name
    printf, unit, 'observation_date ', observation_date
    printf, unit, 'Temperature', ' ', 'Temperature Response '

    for j = 0, n_elements(temp_resp[index[i]].temp_resp) - 1 do begin ; -1
      printf, unit, temp_resp[index[i]].temp[j], ' ', temp_resp[index[i]].temp_resp[j]
    endfor

    close, unit, /force
  endfor

end
