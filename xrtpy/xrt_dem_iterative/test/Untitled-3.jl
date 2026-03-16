
obs_index = ['Be-med', 'Be-thin', 'Al-poly', 'Al-poly/Ti-poly', 'Ti-poly', 'Al-thick']
obs_val   = [9052.93, 37795.0, 136018.0, 53763.7, 66925.8, 320.634]
index = [5, 4, 1, 10, 3, 7]
observation_date = '2007-07-10 13:10:51'
-----------------------------
obs_index = [‘Ti-poly’,’Al-thick’, 'Al-poly’, ‘C-poly’, ‘Be-thin’, ‘Be-med’, ‘Al-med’ ]
  obs_val   = [16696880.61, 597.569475,96577064.12,96577064.12, 4185140.293,17168459.9,8075513.799]  
  index     = [3,7,1,2,4,5,6]                     
  observation_date = '2012-10-27 16:27:46'


-----------------------------
lets shift the focus and do something for simple first- then we'll built the code to do something for complex.
  First, lets just get this code to run as a whole and create the sav file without the outputs and plots as well. 
  Let the first plot be the first DEM and let the 2nd plot be the  MC of all DEMs. 


obs_index = ["Al-mesh", "Ti-poly", "Al-poly", "Be-thin"]
index     = [0,3,1,4]  
obs_val = np.array([799.2547, 119.7911, 699.8442, 14.5166])
observation_date = "2009-07-30T07:40:49"
monte_carlo_runs = 10

print,observation_date,index,obs_index,obs_val

wave_resp = make_xrt_wave_resp(contam_time=observation_date, chn_filename='xrt_channels_v0017')
all_temp_resp = make_xrt_temp_resp(wave_resp, /chianti_default)
temp_resp     = all_temp_resp[index]
xrt_dem_iterative2, obs_index, obs_val, temp_resp, logT, dem


plot, logT, alog10(dem[*,0]), psym=10,title='DEM: ' + observation_date, xtitle='log T (K)', ytitle='log DEM'

xrt_dem_iterative2, obs_index, obs_val, temp_resp, logT_mc, dem_mc, mc_iter = monte_carlo_runs
for i=1,monte_carlo_runs do oplot,logt,alog10(dem[*,i]),psym=10,linestyle=1


window, 1, xsize=800, ysize=600
plot, logT_mc, alog10(dem_mc[*,0]), psym=10, title='Monte Carlo DEM (10 iterations)', xtitle='log T (K)', ytitle='log DEM'

for i=1,monte_carlo_runs do oplot,logt,alog10(dem[*,i]),psym=10,linestyle=1

formatted_date = strreplace(observation_date, ' ', '_')
formatted_date = strreplace(formatted_date, ':', '-')
save_filename = 'xrt_dem_output_' + formatted_date + '_MCITER10.sav'
save, logT_mc, dem_mc, filename=save_filename
print, 'Saved DEM file: ', save_filename
