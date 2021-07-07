import pytest
from astropy import units as u
import numpy as np
from xrtpy.response.channel import Channel
import pkg_resources
import sunpy
import sunpy.map
from sunpy.data import manager
import scipy.io
import sunpy.io.special

channel_names = [
    "Al-mesh",
    "Al-poly",
    "C-poly",
    "Ti-poly",
    "Be-thin",
    "Be-med",
    "Al-med",
    "Al-thick",
    "Be-thick",
    "Al-poly/Al-mesh",
    "Al-poly/Ti-poly",
    "Al-poly/Al-thick",
    "Al-poly/Be-thick",
    "C-poly/Ti-poly",
]


@pytest.mark.parametrize("channel_name", channel_names)
def test_channel_name(channel_name):
    channel = Channel(channel_name)
    assert channel.name == channel_name

filename = pkg_resources.resource_filename( "xrtpy", "data/channels/xrt_channels_v0016.genx")

v6_genx = sunpy.io.special.genx.read_genx(filename)
v6_genx_s = v6_genx['SAVEGEN0']

_channel_name_to_index_mapping = {
    "Al-mesh": 0, 
    "Al-poly": 1, 
    "C-poly": 2 ,
    "Ti-poly": 3,
    "Be-thin": 4,
    "Be-med": 5,
    "Al-med": 6,
    "Al-thick": 7, 
    "Be-thick": 8 ,
    "Al-poly/Al-mesh": 9,
    "Al-poly/Ti-poly": 10, 
    "Al-poly/Al-thick": 11,
    "Al-poly/Be-thick": 12 , 
    "C-poly/Ti-poly": 13, 
}


def test_CCD_units():
    wavelength_unit_name = v6_genx_s[0]['CCD']['WAVE_UNITS']
    if wavelength_unit_name == 'Angstroms':
        assert wavelength_unit_name == 'Angstroms'
        wavelength_unit_name = 'Angstrom'
        wavelength_CCD_unit = u.Unit(wavelength_unit_name)
        return wavelength_CCD_unit

    pixel_size_unit_name_IDL = v6_genx_s[0]['CCD']['PIXEL_SIZE_UNITS']
    pixel_size_unit_name = u.Unit(pixel_size_unit_name_IDL)

@pytest.mark.parametrize("channel_name", channel_names)
def test_CCD_wavelength(channel_name):
    #for i in _channel_name_to_index_mapping:
    channel_filter =  Channel(channel_name)

    ccd_wavelength_length = int(channel_filter.ccd.number_of_wavelengths)   #########Repeated######
    ccd_wavelength = channel_filter.ccd.ccd_wavelength[:ccd_wavelength_length]

    IDL_ccd_array_length = int(v6_genx_s[_channel_name_to_index_mapping[channel_name]]['CCD']["LENGTH"])  #########Repeated######
    IDL_ccd_wavelength_AUTO = v6_genx_s[_channel_name_to_index_mapping[channel_name]]['CCD']["WAVE"][:IDL_ccd_array_length] * u.angstrom
    
    assert u.allclose(IDL_ccd_wavelength_AUTO, ccd_wavelength)

    # All filters have the same starting values 
    IDL_ccd_wavelength_MANU = [1.00000,1.10000,1.20000,1.30000,1.40000, 1.50000,1.60000,1.70000,1.80000,1.90000] * u.angstrom
    assert u.allclose(IDL_ccd_wavelength_MANU,ccd_wavelength[0:10])

def test_CCD_quantum_efficiency():
    for i in _channel_name_to_index_mapping:
        channel_filter =  Channel(i)

        ccd_array_length = int(channel_filter.ccd.number_of_wavelengths) 
        ccd_quantum_efficiency = channel_filter.ccd.ccd_quantum_efficiency[:ccd_array_length]
        
        #IDL_ccd_quantum_efficiency_AUTO = v6_genx_s[1]['CCD']["QE"]
        IDL_ccd_array_length = int(v6_genx_s[1]['CCD']["LENGTH"])
        IDL_ccd_quantum_efficiency_AUTO = v6_genx_s[1]['CCD']["QE"][:IDL_ccd_array_length]
            
        assert u.allclose(IDL_ccd_quantum_efficiency_AUTO, ccd_quantum_efficiency)

def test_CCD_pixel_size():
    for i in _channel_name_to_index_mapping:
        
        channel_filter =  Channel(i)
        ccd_pixel_size = channel_filter.ccd.ccd_pixel_size 

        IDL_ccd_quantum_efficiency_AUTO = v6_genx_s[ _channel_name_to_index_mapping[i] ]['CCD']["PIXEL_SIZE"]* u.micron

        assert u.allclose(IDL_ccd_quantum_efficiency_AUTO,ccd_pixel_size)

        
def test_ccd_gain_left():
    for i in _channel_name_to_index_mapping:
        channel_filter =  Channel(i)
        ccd_gain_left = channel_filter.ccd.ccd_gain_left

        IDL_ccd_gain_left_AUTO = v6_genx_s[ _channel_name_to_index_mapping[i] ]['CCD']["GAIN_L"]* u.electron

        assert u.isclose(ccd_gain_left,IDL_ccd_gain_left_AUTO)

def test_ccd_gain_right():
    for i in _channel_name_to_index_mapping:
        channel_filter =  Channel(i)
        ccd_gain_right = channel_filter.ccd.ccd_gain_right

        IDL_ccd_gain_right_AUTO = v6_genx_s[ _channel_name_to_index_mapping[i] ]['CCD']["GAIN_R"]* u.electron

        if ccd_gain_right == IDL_ccd_gain_right_AUTO:
            pass
        else:
            raise ValueError("FAIL: test_ccd_gain_right")

def test_ccd_full_well():
    for i in _channel_name_to_index_mapping:
        channel_filter =  Channel(i)
        ccd_full_well = channel_filter.ccd.ccd_full_well

        IDL_ccd_full_well_AUTO = v6_genx_s[ _channel_name_to_index_mapping[i] ]['CCD']["FULL_WELL"]* u.electron
        
        if ccd_full_well == IDL_ccd_full_well_AUTO:
            pass
        else:
            raise ValueError("FAIL: test_ccd_full_well")

def test_ccd_ev_ore_electron():
    for i in _channel_name_to_index_mapping:
        channel_filter =  Channel(i)
        ccd_full_well = channel_filter.ccd.ccd_ev_ore_electron
 
        IDL_ccd_full_well_AUTO = v6_genx_s[ _channel_name_to_index_mapping[i] ]['CCD']["EV_PER_EL"]* (u.eV/u.electron)

        if ccd_full_well == IDL_ccd_full_well_AUTO:
            pass
        else:
            raise ValueError("FAIL: test_ccd_ev_ore_electron")


def test_ccd_name():
    for i in _channel_name_to_index_mapping:
        channel_filter =  Channel(i)
        ccd_name = channel_filter.ccd.ccd_name

        IDL_ccd_name_AUTO = v6_genx_s[ _channel_name_to_index_mapping[i] ]['CCD']["LONG_NAME"]

        if ccd_name == IDL_ccd_name_AUTO:
            pass
        else:
            raise ValueError("FAIL: test_ccd_name")

########################### EN Filter ###########################

def test_entrancefilter_name():
    for i in _channel_name_to_index_mapping:
        channel_filter =  Channel(i)
        entrancefilter_name = channel_filter.entrancefilter.entrancefilter_name

        IDL_entrancefilter_name_AUTO = v6_genx_s[ _channel_name_to_index_mapping[i] ]['EN_FILTER']["LONG_NAME"]

        if entrancefilter_name == IDL_entrancefilter_name_AUTO:
            pass
        else:
            raise ValueError("FAIL: test_entrancefilter_name")
        

def test_entrancefilter_material():
    for i in _channel_name_to_index_mapping:
        channel_filter =  Channel(i)
        entrancefilter_material = channel_filter.entrancefilter.entrancefilter_material
 
        IDL_entrancefilter_material_AUTO = v6_genx_s[ _channel_name_to_index_mapping[i] ]['EN_FILTER']["MATERIAL"]

    if all( entrancefilter_material == IDL_entrancefilter_material_AUTO):
        pass
    else:
        raise ValueError("FAIL: test_entrancefilter_material")

def test_entrancefilter_thickness():
    for i in _channel_name_to_index_mapping:
        channel_filter =  Channel(i)
        entrancefilter_thickness = channel_filter.entrancefilter.entrancefilter_thickness

        IDL_entrancefilter_thick_AUTO = v6_genx_s[ _channel_name_to_index_mapping[i] ]['EN_FILTER']["THICK"]*u.angstrom

        if np.all(entrancefilter_thickness == IDL_entrancefilter_thick_AUTO):
            pass
        else:
            raise ValueError("FAIL: test_entrancefilter_thickness")

def test_entrancefilter_density():
    for i in _channel_name_to_index_mapping:
        channel_filter =  Channel(i)
        entrancefilter_density = channel_filter.entrancefilter.entrancefilter_density

        IDL_entrancefilter_density_AUTO = v6_genx_s[ _channel_name_to_index_mapping[i] ]['EN_FILTER']["DENS"]*(u.g * u.cm**-3)
        
        if np.all(entrancefilter_density == IDL_entrancefilter_density_AUTO):
            pass
        else:
            raise ValueError("FAIL: test_entrancefilter_density")

def test_entrancefilter_wavelength():
    for i in _channel_name_to_index_mapping:
        channel_filter =  Channel(i)

        entrancefilter_wavelength_length = int(channel_filter.entrancefilter.number_of_wavelengths)  #########Repeated######
        entrancefilter_wavelength = channel_filter.entrancefilter.entrancefilter_wavelength[:entrancefilter_wavelength_length]

        IDL_entrancefilter_array_length = int(v6_genx_s[_channel_name_to_index_mapping[i]]['EN_FILTER']["LENGTH"])  #########Repeated######
        IDL_entrancefilter_wavelength_AUTO = v6_genx_s[_channel_name_to_index_mapping[i]]['EN_FILTER']["WAVE"][:IDL_entrancefilter_array_length] * u.Unit('Angstrom')#wavelength_CCD_unit

        assert u.allclose(IDL_entrancefilter_wavelength_AUTO, entrancefilter_wavelength)

def test_entrancefilter_transmission():
    for i in _channel_name_to_index_mapping:
        channel_filter =  Channel(i)

        entrancefilter_transmission_length = int(channel_filter.entrancefilter.number_of_wavelengths)  #########Repeated######
        entrancefilter_transmission = channel_filter.entrancefilter.entrancefilter_transmission[:entrancefilter_transmission_length]

        IDL_entrancefilter_array_length = int(v6_genx_s[_channel_name_to_index_mapping[i]]['EN_FILTER']["LENGTH"])  #########Repeated######
        IDL_entrancefilter_transmission_AUTO = v6_genx_s[_channel_name_to_index_mapping[i]]['EN_FILTER']["TRANS"][:IDL_entrancefilter_array_length] 

        assert u.allclose(IDL_entrancefilter_transmission_AUTO, entrancefilter_transmission)

def test_entrancefilter_mesh_transmission():
    for i in _channel_name_to_index_mapping:
        channel_filter =  Channel(i)
        entrancefilter_mesh_transmission = channel_filter.entrancefilter.entrancefilter_mesh_transmission

        IDL_entrancefilter_mesh_transmission_AUTO = v6_genx_s[ _channel_name_to_index_mapping[i] ]['EN_FILTER']["MESH_TRANS"]
        
        if entrancefilter_mesh_transmission == IDL_entrancefilter_mesh_transmission_AUTO:
            pass
        else:
            raise ValueError("FAIL: test_entrancefilter_mesh_transmission")

def test_entrancefilter_substrate():
    for i in _channel_name_to_index_mapping:
        channel_filter =  Channel(i)
        entrancefilter_substrate = channel_filter.entrancefilter.entrancefilter_substrate

        IDL_entrancefilter_substrate_AUTO = v6_genx_s[ _channel_name_to_index_mapping[i] ]['EN_FILTER']["SUBSTRATE"]
        
        if entrancefilter_substrate == IDL_entrancefilter_substrate_AUTO:
            pass
        else:
            raise ValueError("FAIL: test_entrancefilter_substrate")



########################### Filter_1 ###########################

def test_filter1_name():
    for i in _channel_name_to_index_mapping:
        channel_filter =  Channel(i)
        filter_name = channel_filter.filter_1.filter_name

        IDL_filter_name_AUTO = v6_genx_s[ _channel_name_to_index_mapping[i] ]['FP_FILTER1']["LONG_NAME"]

        if filter_name == IDL_filter_name_AUTO:
            pass
        else:
            raise ValueError("FAIL: test_filter1_name")

def test_filter1_material():
    for i in _channel_name_to_index_mapping:
        channel_filter =  Channel(i)
        filter_material = channel_filter.filter_1.filter_material
 
        IDL_filter_material_AUTO = v6_genx_s[ _channel_name_to_index_mapping[i] ]['FP_FILTER1']["MATERIAL"]

    if all( filter_material == IDL_filter_material_AUTO):
        pass
    else:
        raise ValueError("FAIL: test_filter1_material")

def test_filter1_thickness():
    for i in _channel_name_to_index_mapping:
        channel_filter =  Channel(i)
        filter_thickness = channel_filter.filter_1.filter_thickness

        IDL_filter_thick_AUTO = v6_genx_s[ _channel_name_to_index_mapping[i] ]['FP_FILTER1']["THICK"]*u.angstrom

        if np.all(filter_thickness == IDL_filter_thick_AUTO):
            pass
        else:
            raise ValueError("FAIL: test_filter1_thickness")

def test_filter1_density():
    for i in _channel_name_to_index_mapping:
        channel_filter =  Channel(i)
        filter_density = channel_filter.filter_1.filter_density

        IDL_filter_density_AUTO = v6_genx_s[ _channel_name_to_index_mapping[i] ]['FP_FILTER1']["DENS"]*(u.g * u.cm**-3)
        
        if np.all(filter_density == IDL_filter_density_AUTO):
            pass
        else:
            raise ValueError("FAIL: test_filter1_density")

def test_filter1_wavelength():
    for i in _channel_name_to_index_mapping:
        channel_filter =  Channel(i)

        filter_wavelength_length = int(channel_filter.filter_1.number_of_wavelengths)  #########Repeated######
        filter_wavelength = channel_filter.filter_1.filter_wavelength[:filter_wavelength_length]

        IDL_filter_array_length = int(v6_genx_s[_channel_name_to_index_mapping[i]]['FP_FILTER1']["LENGTH"])  #########Repeated######
        IDL_filter_wavelength_AUTO = v6_genx_s[_channel_name_to_index_mapping[i]]['FP_FILTER1']["WAVE"][:IDL_filter_array_length] * u.Unit('Angstrom')#wavelength_CCD_unit

        assert u.allclose(IDL_filter_wavelength_AUTO, filter_wavelength)

def test_filter1_transmission():
    for i in _channel_name_to_index_mapping:
        channel_filter =  Channel(i)

        filter_transmission_length = int(channel_filter.filter_1.number_of_wavelengths)  #########Repeated######
        filter_transmission = channel_filter.filter_1.filter_transmission[:filter_transmission_length]

        IDL_filter_array_length = int(v6_genx_s[_channel_name_to_index_mapping[i]]['FP_FILTER1']["LENGTH"])  #########Repeated######
        IDL_filter_transmission_AUTO = v6_genx_s[_channel_name_to_index_mapping[i]]['FP_FILTER1']["TRANS"][:IDL_filter_array_length] 

        assert u.allclose(IDL_filter_transmission_AUTO, filter_transmission)

def test_filter1_mesh_transmission():
    for i in _channel_name_to_index_mapping:
        channel_filter =  Channel(i)
        filter_mesh_transmission = channel_filter._filter_1.filter_mesh_trans

        IDL_filter_mesh_transmission_AUTO = v6_genx_s[ _channel_name_to_index_mapping[i] ]['FP_FILTER1']["MESH_TRANS"]
        
        if filter_mesh_transmission == IDL_filter_mesh_transmission_AUTO:
            pass
        else:
            raise ValueError("FAIL: test_filter1_mesh_transmission")

def test_filter1_substrate():
    for i in _channel_name_to_index_mapping:
        channel_filter =  Channel(i)
        filter_substrate = channel_filter.filter_1.filter_substrate

        IDL_filter_substrate_AUTO = v6_genx_s[ _channel_name_to_index_mapping[i] ]['FP_FILTER1']["SUBSTRATE"]
        
        if filter_substrate == IDL_filter_substrate_AUTO:
            pass
        else:
            raise ValueError("FAIL: test_filter1_substrate")

########################### Filter_2 ###########################
def test_filter2_name():
    for i in _channel_name_to_index_mapping:
        channel_filter =  Channel(i)
        filter_name = channel_filter.filter_2.filter_name

        IDL_filter_name_AUTO = v6_genx_s[ _channel_name_to_index_mapping[i] ]['FP_FILTER2']["LONG_NAME"]

        if filter_name == IDL_filter_name_AUTO:
            pass
        else:
            raise ValueError("FAIL: test_filter2_name")

def test_filter2_material():
    for i in _channel_name_to_index_mapping:
        channel_filter =  Channel(i)
        filter_material = channel_filter.filter_2.filter_material
 
        IDL_filter_material_AUTO = v6_genx_s[ _channel_name_to_index_mapping[i] ]['FP_FILTER2']["MATERIAL"]

    if all( filter_material == IDL_filter_material_AUTO):
        pass
    else:
        raise ValueError("FAIL: test_filter2_material")



def test_filter2_thickness():
    for i in _channel_name_to_index_mapping:
        channel_filter =  Channel(i)
        filter_thickness = channel_filter.filter_2.filter_thickness

        IDL_filter_thick_AUTO = v6_genx_s[ _channel_name_to_index_mapping[i] ]['FP_FILTER2']["THICK"]*u.angstrom

        if np.all(filter_thickness == IDL_filter_thick_AUTO):
            pass
        else:
            raise ValueError("FAIL: test_filter2_thickness")

def test_filter2_density():
    for i in _channel_name_to_index_mapping:
        channel_filter =  Channel(i)
        filter_density = channel_filter.filter_2.filter_density

        IDL_filter_density_AUTO = v6_genx_s[ _channel_name_to_index_mapping[i] ]['FP_FILTER2']["DENS"]*(u.g * u.cm**-3)
        
        if np.all(filter_density == IDL_filter_density_AUTO):
            pass
        else:
            raise ValueError("FAIL: test_filter2_density")

def test_filter2_wavelength():
    for i in _channel_name_to_index_mapping:
        channel_filter =  Channel(i)

        filter_wavelength_length = int(channel_filter.filter_2.number_of_wavelengths)  #########Repeated######
        filter_wavelength = channel_filter.filter_2.filter_wavelength[:filter_wavelength_length]

        IDL_filter_array_length = int(v6_genx_s[_channel_name_to_index_mapping[i]]['FP_FILTER2']["LENGTH"])  #########Repeated######
        IDL_filter_wavelength_AUTO = v6_genx_s[_channel_name_to_index_mapping[i]]['FP_FILTER2']["WAVE"][:IDL_filter_array_length] * u.Unit('Angstrom')#wavelength_CCD_unit

        assert u.allclose(IDL_filter_wavelength_AUTO, filter_wavelength)

def test_filter2_transmission():
    for i in _channel_name_to_index_mapping:
        channel_filter =  Channel(i)

        filter_transmission_length = int(channel_filter.filter_2.number_of_wavelengths)  #########Repeated######
        filter_transmission = channel_filter.filter_2.filter_transmission[:filter_transmission_length]

        IDL_filter_array_length = int(v6_genx_s[_channel_name_to_index_mapping[i]]['FP_FILTER2']["LENGTH"])  #########Repeated######
        IDL_filter_transmission_AUTO = v6_genx_s[_channel_name_to_index_mapping[i]]['FP_FILTER2']["TRANS"][:IDL_filter_array_length] 

        assert u.allclose(IDL_filter_transmission_AUTO, filter_transmission)

def test_filter2_mesh_transmission():
    for i in _channel_name_to_index_mapping:
        channel_filter =  Channel(i)
        filter_mesh_transmission = channel_filter.filter_2.filter_mesh_trans

        IDL_filter_mesh_transmission_AUTO = v6_genx_s[ _channel_name_to_index_mapping[i] ]['FP_FILTER2']["MESH_TRANS"]
        
        if filter_mesh_transmission == IDL_filter_mesh_transmission_AUTO:
            pass
        else:
            raise ValueError("FAIL: test_filter2_mesh_transmission")

def test_filter2_substrate():
    for i in _channel_name_to_index_mapping:
        channel_filter =  Channel(i)
        filter_substrate = channel_filter.filter_2.filter_substrate

        IDL_filter_substrate_AUTO = v6_genx_s[ _channel_name_to_index_mapping[i] ]['FP_FILTER2']["SUBSTRATE"]
        
        if filter_substrate == IDL_filter_substrate_AUTO:
            pass
        else:
            raise ValueError("FAIL: test_filter2_substrate")


########################### Geomerty  ###########################

def test_geometry_name():
    for i in _channel_name_to_index_mapping:
        channel_filter =  Channel(i)
        geometry_name = channel_filter.geometry.name

        IDL_geometry_name_AUTO = v6_genx_s[ _channel_name_to_index_mapping[i] ]['GEOM']["LONG_NAME"]

        if geometry_name == IDL_geometry_name_AUTO:
            pass
        else:
            raise ValueError("FAIL: test_geometry_name")



def test_geometry_focal_len():  
    for i in _channel_name_to_index_mapping:
        channel_filter =  Channel(i)
        geometry_focal_len = channel_filter.geometry.focal_len

        IDL_geometry_focal_len_AUTO = v6_genx_s[ _channel_name_to_index_mapping[i] ]['GEOM']["FOC_LEN"]*u.cm
        
        if geometry_focal_len == IDL_geometry_focal_len_AUTO:
            pass
        else:
            raise ValueError("FAIL: test_geometry_focal_len")

def test_geometry_aperture_area():  
    for i in _channel_name_to_index_mapping:
        channel_filter =  Channel(i)
        geometry_aperture_area = channel_filter.geometry.aperture_area

        IDL_geometry_aperture_area_AUTO = v6_genx_s[ _channel_name_to_index_mapping[i] ]['GEOM']["APERTURE_AREA"]*u.cm**2
        
        if geometry_aperture_area == IDL_geometry_aperture_area_AUTO:
            pass
        else:
            raise ValueError("FAIL: test_geometry_aperture_area")

########################### Mirror_1  ###########################

def test_mirror1_name():
    for i in _channel_name_to_index_mapping:
        channel_filter =  Channel(i)
        mirror_name = channel_filter.mirror_1.mirror_name

        IDL_mirror_name_AUTO = v6_genx_s[ _channel_name_to_index_mapping[i] ]['MIRROR1']["LONG_NAME"]

        if mirror_name == IDL_mirror_name_AUTO:
            pass
        else:
            raise ValueError("FAIL: test_mirror1_name")

def test_mirror1_material():
    for i in _channel_name_to_index_mapping:
        channel_filter =  Channel(i)
        mirror_material = channel_filter.mirror_1.mirror_material
 
        IDL_mirror_material_AUTO = v6_genx_s[ _channel_name_to_index_mapping[i] ]['MIRROR1']["MATERIAL"]

    if mirror_material == IDL_mirror_material_AUTO:
        pass
    else:
        raise ValueError("FAIL: test_mirror1_material")

def test_mirror1_density():
    for i in _channel_name_to_index_mapping:
        channel_filter =  Channel(i)
        mirror_density = channel_filter.mirror_1.mirror_density

        IDL_mirror_density_AUTO = v6_genx_s[ _channel_name_to_index_mapping[i] ]['MIRROR1']["DENS"]*(u.g * u.cm**-3)
        
        if np.all(mirror_density == IDL_mirror_density_AUTO):
            pass
        else:
            raise ValueError("FAIL: test_mirror1_density")

def test_mirro1_graze_angle():
    for i in _channel_name_to_index_mapping:
        channel_filter =  Channel(i)
        mirror_graze_angle = channel_filter.mirror_1.mirror_graze_angle

        IDL_mirror_graze_angle_AUTO = v6_genx_s[ _channel_name_to_index_mapping[i] ]['MIRROR1']['GRAZE_ANGLE']* u.deg
        
        if mirror_graze_angle == IDL_mirror_graze_angle_AUTO:
            pass
        else:
            raise ValueError("FAIL: test_mirro1_graze_angle")


def test_mirror1_wavelength():
    for i in _channel_name_to_index_mapping:
        channel_filter =  Channel(i)

        mirror_number_of_length = int(channel_filter.mirror_1.number_of_wavelengths)
        mirror_wavelength = channel_filter.mirror_1.mirror_wavelength[:mirror_number_of_length]

        IDL_mirror_array_length = int(v6_genx_s[_channel_name_to_index_mapping[i]]['MIRROR1']["LENGTH"])  #########Repeated######
        IDL_mirror_wavelength_AUTO = v6_genx_s[_channel_name_to_index_mapping[i]]['MIRROR1']["WAVE"][:IDL_mirror_array_length] * u.Unit('Angstrom')#wavelength_CCD_unit

        assert u.allclose(IDL_mirror_wavelength_AUTO, mirror_wavelength)


def test_mirror1_reflection():
    for i in _channel_name_to_index_mapping:
        channel_filter =  Channel(i)

        mirror_number_of_length = int(channel_filter.mirror_1.number_of_wavelengths)
        mirror_reflection1 = channel_filter.mirror_1.mirror_reflection1[:mirror_number_of_length]

        IDL_mirror_array_length = int(v6_genx_s[_channel_name_to_index_mapping[i]]['MIRROR1']["LENGTH"])  #########Repeated######
        IDL_mirror_wavelength_AUTO = v6_genx_s[_channel_name_to_index_mapping[i]]['MIRROR1']["REFL"][:IDL_mirror_array_length] * u.Unit('Angstrom')#wavelength_CCD_unit

        assert u.allclose(IDL_mirror_wavelength_AUTO, mirror_reflection1)


########################### Mirror_2  ###########################

def test_mirror2_name():
    for i in _channel_name_to_index_mapping:
        channel_filter =  Channel(i)
        mirror_name = channel_filter.mirror_1.mirror_name

        IDL_mirror_name_AUTO = v6_genx_s[ _channel_name_to_index_mapping[i] ]['MIRROR2']["LONG_NAME"]

        if mirror_name == IDL_mirror_name_AUTO:
            pass
        else:
            raise ValueError("FAIL: test_mirror2_name")

def test_mirror2_material():
    for i in _channel_name_to_index_mapping:
        channel_filter =  Channel(i)
        mirror_material = channel_filter.mirror_1.mirror_material
 
        IDL_mirror_material_AUTO = v6_genx_s[ _channel_name_to_index_mapping[i] ]['MIRROR2']["MATERIAL"]

    if mirror_material == IDL_mirror_material_AUTO:
        pass
    else:
        raise ValueError("FAIL: test_mirror2_material")

def test_mirror2_density():
    for i in _channel_name_to_index_mapping:
        channel_filter =  Channel(i)
        mirror_density = channel_filter.mirror_1.mirror_density

        IDL_mirror_density_AUTO = v6_genx_s[ _channel_name_to_index_mapping[i] ]['MIRROR2']["DENS"]*(u.g * u.cm**-3)
        
        if np.all(mirror_density == IDL_mirror_density_AUTO):
            pass
        else:
            raise ValueError("FAIL: test_mirror2_density")

def test_mirror2_graze_angle():
    for i in _channel_name_to_index_mapping:
        channel_filter =  Channel(i)
        mirror_graze_angle = channel_filter.mirror_1.mirror_graze_angle

        IDL_mirror_graze_angle_AUTO = v6_genx_s[ _channel_name_to_index_mapping[i] ]['MIRROR2']['GRAZE_ANGLE']* u.deg
        
        if mirror_graze_angle == IDL_mirror_graze_angle_AUTO:
            pass
        else:
            raise ValueError("FAIL: test_mirror2_graze_angle")

def test_mirror2_wavelength():
    for i in _channel_name_to_index_mapping:
        channel_filter =  Channel(i)

        mirror_number_of_length = int(channel_filter.mirror_1.number_of_wavelengths)
        mirror_wavelength = channel_filter.mirror_1.mirror_wavelength[:mirror_number_of_length]

        IDL_mirror_array_length = int(v6_genx_s[_channel_name_to_index_mapping[i]]['MIRROR2']["LENGTH"])  #########Repeated######
        IDL_mirror_wavelength_AUTO = v6_genx_s[_channel_name_to_index_mapping[i]]['MIRROR2']["WAVE"][:IDL_mirror_array_length] * u.Unit('Angstrom')#wavelength_CCD_unit

        assert u.allclose(IDL_mirror_wavelength_AUTO, mirror_wavelength)


def test_mirror2_reflection():
    for i in _channel_name_to_index_mapping:
        channel_filter =  Channel(i)

        mirror_number_of_length = int(channel_filter.mirror_1.number_of_wavelengths)
        mirror_reflection1 = channel_filter.mirror_1.mirror_reflection1[:mirror_number_of_length]

        IDL_mirror_array_length = int(v6_genx_s[_channel_name_to_index_mapping[i]]['MIRROR2']["LENGTH"])  #########Repeated######
        IDL_mirror_wavelength_AUTO = v6_genx_s[_channel_name_to_index_mapping[i]]['MIRROR2']["REFL"][:IDL_mirror_array_length] * u.Unit('Angstrom')#wavelength_CCD_unit

        assert u.allclose(IDL_mirror_wavelength_AUTO, mirror_reflection1)

########################### Name ###########################

def test_channel_name():
    for i in _channel_name_to_index_mapping:
        channel_filter =  Channel(i)
        name = channel_filter.name

        IDL_mirror_name_AUTO = v6_genx_s[ _channel_name_to_index_mapping[i] ]['NAME']

        if name == IDL_mirror_name_AUTO:
            pass
        else:
            raise ValueError("FAIL: test_channel_name")


########################### ###########################

def test_channel_wavelength():
    for i in _channel_name_to_index_mapping:
        channel_filter =  Channel(i)

        wavelength_length = int(channel_filter.number_of_wavelengths)   #########Repeated######
        wavelength = channel_filter.wavelength[:wavelength_length]

        IDL_array_length = int(v6_genx_s[_channel_name_to_index_mapping[i]]["LENGTH"])  #########Repeated######
        IDL_wavelength_AUTO = v6_genx_s[_channel_name_to_index_mapping[i]]["WAVE"][:IDL_array_length] * u.Unit('Angstrom')#wavelength_CCD_unit
        
        assert u.allclose(IDL_wavelength_AUTO, wavelength)

########################### ###########################

def test_channel_transmission():
    for i in _channel_name_to_index_mapping:
        channel_filter =  Channel(i)

        transmission_length = int(channel_filter.number_of_wavelengths)  #########Repeated######
        transmission = channel_filter.transmission[:transmission_length]

        IDL_array_length = int(v6_genx_s[_channel_name_to_index_mapping[i]]["LENGTH"])  #########Repeated######
        IDL_transmission_AUTO = v6_genx_s[_channel_name_to_index_mapping[i]]["TRANS"][:IDL_array_length] 

        assert u.allclose(IDL_transmission_AUTO, transmission)
        
########################### ###########################

def test_channel_number_of_wavelengths():
    for i in _channel_name_to_index_mapping:
        channel_filter =  Channel(i)

        transmission = channel_filter.number_of_wavelengths

        IDL_array_length = v6_genx_s[_channel_name_to_index_mapping[i]]["LENGTH"]
    if transmission == IDL_array_length:
        pass
    else:
        raise ValueError("FAIL: test_channel_number_of_wavelengths")

########################### ###########################       

def test_channel_observatory():
    for i in _channel_name_to_index_mapping:
        channel_filter =  Channel(i)
        
        observatory = channel_filter.observatory
        
        IDL_observatory = v6_genx_s[_channel_name_to_index_mapping[i]]["OBSERVATORY"]
    
    if observatory == IDL_observatory:
        pass
    else:
        raise ValueError("FAIL: test_channel_observatory")

    
########################### ###########################       

def test_channel_instrument():
    for i in _channel_name_to_index_mapping:
        channel_filter =  Channel(i)
        
        instrument = channel_filter.instrument
        
        IDL_instrument = v6_genx_s[_channel_name_to_index_mapping[i]]["INSTRUMENT"]
    
    if instrument == IDL_instrument:
        pass
    else:
        raise ValueError("FAIL: test_channel_instrument")
