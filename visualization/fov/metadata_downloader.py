from astropy.io import fits
import astropy.units as u

import numpy as np

from astropy.coordinates import SkyCoord

from sunpy.coordinates import frames
from astropy.time import Time
import scipy.io as sio
from scipy.io import readsav
import urllib.request


def download_metadata(xrt_downloaded_files, filen, overwrite=False):
    url_lis = xrt_downloaded_files[0][:]['fileid']
    url_str_lis = []
    for i in range(len(url_lis)):
        url_str_lis.append(url_lis[i])
    primary_hdu = fits.PrimaryHDU(data=np.ones((3, 3)))
    c1 = fits.Column(name='URL', array=url_str_lis, format='100A')
    c2 = fits.Column(name='header_int', array=np.asarray(range(len(url_str_lis)))+2, format='J')
    table_hdu = fits.BinTableHDU.from_columns([c1, c2])

    hdul2 = fits.HDUList([primary_hdu, table_hdu])
    #hlis = []

    for i in range(len(url_str_lis)):
        if (i%10 == 0):
            print(int(1000.0*i/len(url_str_lis))/10.0,'%')
        fsspec_kwargs = {"block_size": 100_000, "cache_type": "bytes"}
        with fits.open(url_lis[i], use_fsspec=True, fsspec_kwargs=fsspec_kwargs) as hdul:
            #hlis.append(hdul[0].header)
            # Download a single header
            t_header = hdul[0].header
            image_hdu = fits.ImageHDU(data=np.ones((100, 100)), header=t_header,name="header"+str(i))
        hdul2.append(image_hdu)        
    return hdul2

    #hdul2.writeto(filen,overwrite=overwrite)

def date_to_meta(xrt_download_list):
    time_lis = xrt_download_list[0]['Start Time']
    year_lis = xrt_download_list[0]['Start Time'].ymdhms.year
    month_lis = xrt_download_list[0]['Start Time'].ymdhms.month
    day_lis = xrt_download_list[0]['Start Time'].ymdhms.day
    new_date = []
    file_lis = []
    for i in range(len(time_lis)):
        year_str  = str(year_lis[i]) 
        month_str = str(month_lis[i])
        day_str   = str(day_lis[i])  
        if len(day_str) < 2:
            day_str = '0'+day_str
        if len(month_str) < 2:
            month_str = '0'+month_str
        ndatei =  year_str + month_str+ day_str  
        if ndatei in file_lis:
            new_date.append(file_lis.index(ndatei))
            #new_date[file_lis.index(ndatei)].append(ndatei)
        else:
            file_lis.append(ndatei)
            new_date.append(file_lis.index(ndatei))
            #print(ndatei,file_lis.index(ndatei))
    return file_lis, new_date

def get_urls(file_n_lis, ggg):
    nfile = len(file_n_lis)
    geny_lis = []
    for i in range(nfile):
        find_url = 'xrt'+file_n_lis[i]
        findex = ggg.find(find_url)
        gen_fn = ggg[findex:findex+35]
        findex2 = gen_fn.find('geny')
        gen_fn = gen_fn[:findex2+4]
        geny_lis.append(gen_fn)
    return geny_lis

def get_metafile(geny_lis):
    url_start = 'https://sot.lmsal.com/data/sot/metadata/sswdb/hinode/xrt/xrt_genxcat/'
    ngeny = len(geny_lis)
    meta_lis = []
    for i in range(ngeny):
        print(i)
        gen_fn = geny_lis[i]
        f, h = urllib.request.urlretrieve(url_start + gen_fn)
        print(i)
        data2 = readsav(f)["p0"]
        data_dict2 = {k : data2[k] for k in data2.dtype.names}
        meta_lis.append(data_dict2)
    return meta_lis

#def mk_meta_header(meta_lis):
        #print(data_dict2['DATE_OBS'])

def meta_to_dict(data_dict, di):
    dkeys = data_dict.keys()
    hdict = {}
    for dki in dkeys:
        try:
            hdict[dki] = data_dict[dki][di].decode('ascii')
        except:
            hdict[dki] = data_dict[dki][di]
    return hdict

def match_vso_to_cat(data_dict_lis, cat_fi, xrt_download):
    
    n_dict = len(data_dict_lis)
    cat_time_lis = []
    for i in range(n_dict):
        data_dict = data_dict_lis[i]
        date_obs_cat = data_dict['DATE_OBS']
        cat_len = len(date_obs_cat)
        cat_str = []
        for cat_bin in date_obs_cat:
            cat_str.append(cat_bin.decode('ascii'))
        cat_time = Time(np.asarray(cat_str), format='isot', scale='utc')
        cat_time_lis.append(cat_time)
    #print('yo')
    min_ti_lis = []
    delt_lis = []
    delt_lisp = []
    delt_lism = []
    header_lis = []
    for i in range(len(xrt_download[0]['Start Time'])):
        cat_time = cat_time_lis[cat_fi[i]]
        stime = xrt_download[0]['Start Time'][i]
        delt = cat_time - stime
        delt = delt.value*24.0*3600.0
        min_ti = np.argmin(np.abs(delt))
        min_ti_lis.append(min_ti)
        delt_lis.append(delt[min_ti])
        try:
            delt_lisp.append(delt[min_ti+1])
            delt_lism.append(delt[min_ti-1])
        except:
            print()
        header_lis.append(meta_to_dict(data_dict_lis[cat_fi[i]],min_ti))
    return header_lis

def get_html_lis():
    url_start = 'https://sot.lmsal.com/data/sot/metadata/sswdb/hinode/xrt/xrt_genxcat/'
    with urllib.request.urlopen(url_start) as response:  
        
        html = response.read()
        #print(response.info())
        #print()
        ggg = html.decode('utf-8')
    return ggg

def download_metadata_fast(xrt_downloaded_files, ggg=None):
    if (ggg == None):
        ggg = get_html_lis()
    file_lis, new_date = date_to_meta(xrt_downloaded_files)
    genyl = get_urls(file_lis, ggg)
    print('downloading')
    tmeta_lis = get_metafile(genyl)
    print('done')
    hlis3 = match_vso_to_cat(tmeta_lis, new_date, xrt_downloaded_files)
    return hlis3

def fetch_metadata(xrt_downloaded_files, fast_bool = True):
    if fast_bool:
        print('Fast Metadata (Level 0)')
        return download_metadata_fast(xrt_downloaded_files, ggg=None)
    else:
        print('Slow Metadata (Level 1)')
        hdul = download_metadata(xrt_downloaded_files,'')
        hlis = []
        for i in range(len(xrt_downloaded_files[0])):
            hlis.append(hdul[i+2].header)
        return hlis
    

