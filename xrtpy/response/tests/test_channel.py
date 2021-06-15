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

@pytest.mark.parametrize("cls, args, attribute", cls_args_attribute)
def test_attributes(cls, args, attribute):
    instance = cls(*args)
    getattr(instance, attribute)

class ExampleClass():
    
    @property
    def first_property(self):
        return 1
    
    @property
    def second_property(self):
        return "2"

#put the names of properties we want to check

properties_to_check = ["first_property", "second_property", "third_property"]

@pytest.mark.parametrize("property", properties_to_check)
def test_accessing_property(property):
    instance = ExampleClass()  # change to name of class being checked
    getattr(instance, property)

# put the names of properties and the expected values
properties_and_expected_values = [
    ("first_property", 1), 
    ("second_property", "2"), 
    ("third_property", "three"), # this test should fail
]

@pytest.mark.parametrize("property, expected_value", properties_and_expected_values)
def test_accessing_property(property, expected_value):
    instance = ExampleClass()  # change to name of class being checked
    actual_value = getattr(instance, property)
    assert actual_value == expected_value