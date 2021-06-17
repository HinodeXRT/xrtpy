import pytest
from xrtpy.response.channel import Channel

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

cls_args_attribute = [
    (
        Channel,  # class
        [
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
        ],  # args
        [
            "name", # first attribute
            "wavelength",  # second attribute
            "transmission",
            "number_of_wavelengths",
            "observatory",
            "instrument",
            "geometry",
            "entrancefilter",
            "mirror",
            "filter",
            "ccd",
        ],
    ),
    (
        Geometry, 
        [],
        [
            "name",
            "focal_len", 
            "aperture_area",  
    ),
    (
        EntranceFilter,  
        [], 
        [
            "entrancefilter_name",
            "entrancefilter_material", 
            "entrancefilter_thickness",
            "entrancefilter_density",
            "entrancefilter_wavelength",
            "entrancefilter_transmission",
            "number_of_wavelengths",
            "entrancefilter_mesh_transmission",
            "entrancefilter_substrate",
        ]
    ),
    (
        Mirror,  
        [],  
        [
            "mirror_name",
            "mirror_material", 
            "mirror_density",  
            "mirror_graze_angle",
            "mirror_wavelength",
            "mirror_reflection1",
            "number_of_wavelengths",
        ]
    ),
    (
        Filter, 
        [],  
        [
            "filter_name",
            "filter_material", 
            "filter_thickness",  
            "filter_density",
            "filter_wavelength", 
            "filter_transmission", 
            "number_of_wavelengths",
            "filter_mesh_trans",  
            "filter_substrate",  
        ]
    ),
    (
        CCD, 
        [],  
        [
            "ccd_name",
            "ccd_ev_ore_electron",  
            "ccd_full_well", 
            "ccd_gain_left",
            "ccd_gain_right",
            "ccd_quantum_efficiency", 
            "number_of_wavelengths",  
            "ccd_pixel_size",  
            "ccd_wavelength", 
        ]
    )  
]


#cls: "Channel,geometry,entrance_filter,mirror,filter,CCD"

#argC name,wavelength,transmission,number_of_wavelengths, observatory, instrument
#argG name, focal_len, aperture_area
#argE name, filter_material,filter_thickness,filter_density, wavelength,transmission,number_of_wavelengths,mesh_trans,substrate
