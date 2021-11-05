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


filename = pkg_resources.resource_filename(
    "xrtpy", "data/channels/xrt_channels_v0016.genx"
)

v6_genx = sunpy.io.special.genx.read_genx(filename)
v6_genx_s = v6_genx["SAVEGEN0"]

_channel_name_to_index_mapping = {
    "Al-mesh": 0,
    "Al-poly": 1,
    "C-poly": 2,
    "Ti-poly": 3,
    "Be-thin": 4,
    "Be-med": 5,
    "Al-med": 6,
    "Al-thick": 7,
    "Be-thick": 8,
    "Al-poly/Al-mesh": 9,
    "Al-poly/Ti-poly": 10,
    "Al-poly/Al-thick": 11,
    "Al-poly/Be-thick": 12,
    "C-poly/Ti-poly": 13,
}


@pytest.mark.parametrize("channel_name", channel_names)
def test_CCD_wavelength(channel_name):
    channel_filter = Channel(channel_name)

    ccd_wavelength_length = int(channel_filter.ccd.number_of_wavelengths)
    ccd_wavelength = channel_filter.ccd.ccd_wavelength[:ccd_wavelength_length]

    idl_ccd_array_length = int(
        v6_genx_s[_channel_name_to_index_mapping[channel_name]]["CCD"]["LENGTH"]
    )
    idl_ccd_wavelength_auto = (
        v6_genx_s[_channel_name_to_index_mapping[channel_name]]["CCD"]["WAVE"][
            :idl_ccd_array_length
        ]
        * u.angstrom
    )

    assert u.allclose(idl_ccd_wavelength_auto, ccd_wavelength)

    idl_ccd_wavelength_manu = [
        1.00000,
        1.10000,
        1.20000,
        1.30000,
        1.40000,
        1.50000,
        1.60000,
        1.70000,
        1.80000,
        1.90000,
    ] * u.angstrom
    assert u.allclose(idl_ccd_wavelength_manu, ccd_wavelength[0:10])


@pytest.mark.parametrize("channel_name", channel_names)
def test_CCD_quantum_efficiency(channel_name):
    channel_filter = Channel(channel_name)

    ccd_array_length = int(channel_filter.ccd.number_of_wavelengths)
    ccd_quantum_efficiency = channel_filter.ccd.ccd_quantum_efficiency[
        :ccd_array_length
    ]

    idl_ccd_array_length = int(
        v6_genx_s[_channel_name_to_index_mapping[channel_name]]["CCD"]["LENGTH"]
    )
    idl_ccd_quantum_efficiency_auto = v6_genx_s[
        _channel_name_to_index_mapping[channel_name]
    ]["CCD"]["QE"][:idl_ccd_array_length]

    assert u.allclose(idl_ccd_quantum_efficiency_auto, ccd_quantum_efficiency)

    idl_ccd_quantum_efficiency_manu = [
        0.0573069,
        0.0751920,
        0.0960381,
        0.119867,
        0.146638,
        0.176252,
        0.208541,
        0.243277,
        0.280167,
        0.318879,
        0.359036,
        0.400219,
        0.441984,
        0.483898,
    ]
    assert idl_ccd_quantum_efficiency_manu, ccd_quantum_efficiency[0:13]


@pytest.mark.parametrize("channel_name", channel_names)
def test_CCD_pixel_size(channel_name):
    channel_filter = Channel(channel_name)
    ccd_pixel_size = channel_filter.ccd.ccd_pixel_size

    idl_ccd_quantum_efficiency_auto = (
        v6_genx_s[_channel_name_to_index_mapping[channel_name]]["CCD"]["PIXEL_SIZE"]
        * u.micron
    )

    assert u.allclose(idl_ccd_quantum_efficiency_auto, ccd_pixel_size)


@pytest.mark.parametrize("channel_name", channel_names)
def test_ccd_gain_left(channel_name):
    channel_filter = Channel(channel_name)
    ccd_gain_left = channel_filter.ccd.ccd_gain_left

    idl_ccd_gain_left_auto = (
        v6_genx_s[_channel_name_to_index_mapping[channel_name]]["CCD"]["GAIN_L"]
        * u.electron
    )

    assert u.isclose(ccd_gain_left, idl_ccd_gain_left_auto)


@pytest.mark.parametrize("channel_name", channel_names)
def test_ccd_gain_right(channel_name):
    channel_filter = Channel(channel_name)
    ccd_gain_right = channel_filter.ccd.ccd_gain_right

    idl_ccd_gain_right_auto = (
        v6_genx_s[_channel_name_to_index_mapping[channel_name]]["CCD"]["GAIN_R"]
        * u.electron
    )

    assert u.isclose(ccd_gain_right, idl_ccd_gain_right_auto)


@pytest.mark.parametrize("channel_name", channel_names)
def test_ccd_full_well(channel_name):
    channel_filter = Channel(channel_name)
    ccd_full_well = channel_filter.ccd.ccd_full_well

    idl_ccd_full_well_auto = (
        v6_genx_s[_channel_name_to_index_mapping[channel_name]]["CCD"]["FULL_WELL"]
        * u.electron
    )

    assert u.isclose(ccd_full_well, idl_ccd_full_well_auto)


@pytest.mark.parametrize("channel_name", channel_names)
def test_ccd_ev_ore_electron(channel_name):

    channel_filter = Channel(channel_name)
    ccd_full_well = channel_filter.ccd.ccd_ev_ore_electron

    idl_ccd_full_well_auto = v6_genx_s[_channel_name_to_index_mapping[channel_name]][
        "CCD"
    ]["EV_PER_EL"] * (u.eV / u.electron)

    assert u.isclose(ccd_full_well, idl_ccd_full_well_auto)


@pytest.mark.parametrize("channel_name", channel_names)
def test_ccd_name(channel_name):
    channel_filter = Channel(channel_name)
    ccd_name = channel_filter.ccd.ccd_name

    idl_ccd_name_auto = v6_genx_s[_channel_name_to_index_mapping[channel_name]]["CCD"][
        "LONG_NAME"
    ]

    assert ccd_name == idl_ccd_name_auto


@pytest.mark.parametrize("channel_name", channel_names)
def test_entrancefilter_name(channel_name):

    channel_filter = Channel(channel_name)
    entrancefilter_name = channel_filter.entrancefilter.entrancefilter_name

    IDL_entrancefilter_name_AUTO = v6_genx_s[
        _channel_name_to_index_mapping[channel_name]
    ]["EN_FILTER"]["LONG_NAME"]

    assert entrancefilter_name == IDL_entrancefilter_name_AUTO


@pytest.mark.parametrize("channel_name", channel_names)
def test_entrancefilter_material(channel_name):

    channel_filter = Channel(channel_name)
    entrancefilter_material = channel_filter.entrancefilter.entrancefilter_material

    idl_entrancefilter_material_auto = v6_genx_s[
        _channel_name_to_index_mapping[channel_name]
    ]["EN_FILTER"]["MATERIAL"]

    if np.all(entrancefilter_material == idl_entrancefilter_material_auto):
        pass
    else:
        raise ValueError("FAIL: test_entrancefilter_material")


@pytest.mark.parametrize("channel_name", channel_names)
def test_entrancefilter_thickness(channel_name):

    channel_filter = Channel(channel_name)
    entrancefilter_thickness = channel_filter.entrancefilter.entrancefilter_thickness

    idl_entrancefilter_thick_auto = (
        v6_genx_s[_channel_name_to_index_mapping[channel_name]]["EN_FILTER"]["THICK"]
        * u.angstrom
    )

    assert u.allclose(entrancefilter_thickness, idl_entrancefilter_thick_auto)


@pytest.mark.parametrize("channel_name", channel_names)
def test_entrancefilter_density(channel_name):
    channel_filter = Channel(channel_name)
    entrancefilter_density = channel_filter.entrancefilter.entrancefilter_density

    idl_entrancefilter_density_auto = v6_genx_s[
        _channel_name_to_index_mapping[channel_name]
    ]["EN_FILTER"]["DENS"] * (u.g * u.cm ** -3)

    assert u.allclose(entrancefilter_density, idl_entrancefilter_density_auto)


@pytest.mark.parametrize("channel_name", channel_names)
def test_entrancefilter_wavelength(channel_name):
    channel_filter = Channel(channel_name)

    entrancefilter_wavelength_length = int(
        channel_filter.entrancefilter.number_of_wavelengths
    )
    entrancefilter_wavelength = channel_filter.entrancefilter.entrancefilter_wavelength[
        :entrancefilter_wavelength_length
    ]

    idl_entrancefilter_array_length = int(
        v6_genx_s[_channel_name_to_index_mapping[channel_name]]["EN_FILTER"]["LENGTH"]
    )
    idl_entrancefilter_wavelength_auto = v6_genx_s[
        _channel_name_to_index_mapping[channel_name]
    ]["EN_FILTER"]["WAVE"][:idl_entrancefilter_array_length] * u.Unit(
        "Angstrom"
    )  # wavelength_CCD_unit

    assert u.allclose(idl_entrancefilter_wavelength_auto, entrancefilter_wavelength)

    idl_entrancefilter_wavelength_manu = [
        1.00000,
        1.00802,
        1.01610,
        1.02424,
        1.03245,
        1.04073,
        1.04907,
        1.05748,
        1.06595,
        1.07450,
    ] * u.angstrom
    assert u.allclose(
        idl_entrancefilter_wavelength_manu, entrancefilter_wavelength[0:10]
    )


@pytest.mark.parametrize("channel_name", channel_names)
def test_entrancefilter_transmission(channel_name):
    channel_filter = Channel(channel_name)

    entrancefilter_transmission_length = int(
        channel_filter.entrancefilter.number_of_wavelengths
    )
    entrancefilter_transmission = (
        channel_filter.entrancefilter.entrancefilter_transmission[
            :entrancefilter_transmission_length
        ]
    )

    idl_entrancefilter_array_length = int(
        v6_genx_s[_channel_name_to_index_mapping[channel_name]]["EN_FILTER"]["LENGTH"]
    )
    idl_entrancefilter_transmission_auto = v6_genx_s[
        _channel_name_to_index_mapping[channel_name]
    ]["EN_FILTER"]["TRANS"][:idl_entrancefilter_array_length]

    assert u.allclose(idl_entrancefilter_transmission_auto, entrancefilter_transmission)


@pytest.mark.parametrize("channel_name", channel_names)
def test_entrancefilter_mesh_transmission(channel_name):
    channel_filter = Channel(channel_name)
    entrancefilter_mesh_transmission = (
        channel_filter.entrancefilter.entrancefilter_mesh_transmission
    )

    idl_entrancefilter_mesh_transmission_auto = v6_genx_s[
        _channel_name_to_index_mapping[channel_name]
    ]["EN_FILTER"]["MESH_TRANS"]

    assert entrancefilter_mesh_transmission == idl_entrancefilter_mesh_transmission_auto


@pytest.mark.parametrize("channel_name", channel_names)
def test_entrancefilter_substrate(channel_name):
    channel_filter = Channel(channel_name)
    entrancefilter_substrate = channel_filter.entrancefilter.entrancefilter_substrate

    idl_entrancefilter_substrate_auto = v6_genx_s[
        _channel_name_to_index_mapping[channel_name]
    ]["EN_FILTER"]["SUBSTRATE"]

    assert entrancefilter_substrate == idl_entrancefilter_substrate_auto


@pytest.mark.parametrize("channel_name", channel_names)
def test_filter1_name(channel_name):
    channel_filter = Channel(channel_name)
    filter_name = channel_filter.filter_1.name

    idl_filter_name_auto = v6_genx_s[_channel_name_to_index_mapping[channel_name]][
        "FP_FILTER1"
    ]["LONG_NAME"]

    assert filter_name == idl_filter_name_auto


@pytest.mark.parametrize("channel_name", channel_names)
def test_filter1_material(channel_name):
    channel_filter = Channel(channel_name)
    filter_material = channel_filter.filter_1.material

    idl_filter_material_auto = v6_genx_s[_channel_name_to_index_mapping[channel_name]][
        "FP_FILTER1"
    ]["MATERIAL"]

    assert np.all(filter_material == idl_filter_material_auto)


@pytest.mark.parametrize("channel_name", channel_names)
def test_filter1_thickness(channel_name):
    channel_filter = Channel(channel_name)
    filter_thickness = channel_filter.filter_1.thickness

    idl_filter_thick_auto = (
        v6_genx_s[_channel_name_to_index_mapping[channel_name]]["FP_FILTER1"]["THICK"]
        * u.angstrom
    )

    assert np.all(filter_thickness == idl_filter_thick_auto)


@pytest.mark.parametrize("channel_name", channel_names)
def test_filter1_density(channel_name):
    channel_filter = Channel(channel_name)
    filter_density = channel_filter.filter_1.density

    idl_filter_density_auto = v6_genx_s[_channel_name_to_index_mapping[channel_name]][
        "FP_FILTER1"
    ]["DENS"] * (u.g * u.cm ** -3)

    assert u.allclose(filter_density, idl_filter_density_auto)


@pytest.mark.parametrize("channel_name", channel_names)
def test_filter1_wavelength(channel_name):
    channel_filter = Channel(channel_name)

    filter_wavelength_length = int(channel_filter.filter_1.number_of_wavelengths)
    filter_wavelength = channel_filter.filter_1.wavelength[:filter_wavelength_length]

    idl_filter_array_length = int(
        v6_genx_s[_channel_name_to_index_mapping[channel_name]]["FP_FILTER1"]["LENGTH"]
    )
    idl_filter_wavelength_auto = (
        v6_genx_s[_channel_name_to_index_mapping[channel_name]]["FP_FILTER1"]["WAVE"][
            :idl_filter_array_length
        ]
        * u.angstrom
    )

    assert u.allclose(idl_filter_wavelength_auto, filter_wavelength)

    idl_filter_wavelength_manu = [
        1.00000,
        1.00802,
        1.01610,
        1.02424,
        1.03245,
        1.04073,
        1.04907,
        1.05748,
        1.06595,
        1.07450,
    ] * u.angstrom
    assert u.allclose(idl_filter_wavelength_manu, filter_wavelength[0:10])


@pytest.mark.parametrize("channel_name", channel_names)
def test_filter1_transmission(channel_name):
    channel_filter = Channel(channel_name)

    filter_transmission_length = int(channel_filter.filter_1.number_of_wavelengths)
    filter_transmission = channel_filter.filter_1.transmission[
        :filter_transmission_length
    ]

    idl_filter_array_length = int(
        v6_genx_s[_channel_name_to_index_mapping[channel_name]]["FP_FILTER1"]["LENGTH"]
    )
    idl_filter_transmission_auto = v6_genx_s[
        _channel_name_to_index_mapping[channel_name]
    ]["FP_FILTER1"]["TRANS"][:idl_filter_array_length]

    assert u.allclose(idl_filter_transmission_auto, filter_transmission)


@pytest.mark.parametrize("channel_name", channel_names)
def test_filter1_mesh_transmission(channel_name):
    channel_filter = Channel(channel_name)
    filter_mesh_transmission = channel_filter._filter_1.mesh_trans

    idl_filter_mesh_transmission_auto = v6_genx_s[
        _channel_name_to_index_mapping[channel_name]
    ]["FP_FILTER1"]["MESH_TRANS"]

    assert filter_mesh_transmission == idl_filter_mesh_transmission_auto


@pytest.mark.parametrize("channel_name", channel_names)
def test_filter1_substrate(channel_name):
    channel_filter = Channel(channel_name)
    filter_substrate = channel_filter.filter_1.substrate

    idl_filter_substrate_auto = v6_genx_s[_channel_name_to_index_mapping[channel_name]][
        "FP_FILTER1"
    ]["SUBSTRATE"]

    assert filter_substrate == idl_filter_substrate_auto


@pytest.mark.parametrize("channel_name", channel_names)
def test_filter2_name(channel_name):
    channel_filter = Channel(channel_name)
    filter_name = channel_filter.filter_2.name

    IDL_filter_name_AUTO = v6_genx_s[_channel_name_to_index_mapping[channel_name]][
        "FP_FILTER2"
    ]["LONG_NAME"]

    assert filter_name == IDL_filter_name_AUTO


@pytest.mark.parametrize("channel_name", channel_names)
def test_filter2_material(channel_name):
    channel_filter = Channel(channel_name)
    filter_material = channel_filter.filter_2.material

    idl_filter_material_auto = v6_genx_s[_channel_name_to_index_mapping[channel_name]][
        "FP_FILTER2"
    ]["MATERIAL"]

    assert np.all(filter_material == idl_filter_material_auto)


@pytest.mark.parametrize("channel_name", channel_names)
def test_filter2_thickness(channel_name):
    channel_filter = Channel(channel_name)
    filter_thickness = channel_filter.filter_2.thickness

    idl_filter_thick_auto = (
        v6_genx_s[_channel_name_to_index_mapping[channel_name]]["FP_FILTER2"]["THICK"]
        * u.angstrom
    )

    assert u.allclose(filter_thickness, idl_filter_thick_auto)


@pytest.mark.parametrize("channel_name", channel_names)
def test_filter2_density(channel_name):
    channel_filter = Channel(channel_name)
    filter_density = channel_filter.filter_2.density

    IDL_filter_density_AUTO = v6_genx_s[_channel_name_to_index_mapping[channel_name]][
        "FP_FILTER2"
    ]["DENS"] * (u.g * u.cm ** -3)

    np.allclose(filter_density, IDL_filter_density_AUTO)


@pytest.mark.parametrize("channel_name", channel_names)
def test_filter2_wavelength(channel_name):
    channel_filter = Channel(channel_name)

    filter_wavelength_length = int(channel_filter.filter_2.number_of_wavelengths)
    filter_wavelength = channel_filter.filter_2.wavelength[:filter_wavelength_length]

    idl_filter_array_length = int(
        v6_genx_s[_channel_name_to_index_mapping[channel_name]]["FP_FILTER2"]["LENGTH"]
    )
    idl_filter_wavelength_auto = (
        v6_genx_s[_channel_name_to_index_mapping[channel_name]]["FP_FILTER2"]["WAVE"][
            :idl_filter_array_length
        ]
        * u.angstrom
    )

    assert u.allclose(idl_filter_wavelength_auto, filter_wavelength)

    idl_filter_wavelength_manu = [
        1.00000,
        1.00802,
        1.01610,
        1.02424,
        1.03245,
        1.04073,
        1.04907,
        1.05748,
        1.06595,
        1.07450,
    ] * u.angstrom
    assert u.allclose(idl_filter_wavelength_manu, filter_wavelength[0:10])


@pytest.mark.parametrize("channel_name", channel_names)
def test_filter2_transmission(channel_name):

    channel_filter = Channel(channel_name)

    filter_transmission_length = int(channel_filter.filter_2.number_of_wavelengths)
    filter_transmission = channel_filter.filter_2.transmission[
        :filter_transmission_length
    ]

    idl_filter_array_length = int(
        v6_genx_s[_channel_name_to_index_mapping[channel_name]]["FP_FILTER2"]["LENGTH"]
    )
    idl_filter_transmission_auto = v6_genx_s[
        _channel_name_to_index_mapping[channel_name]
    ]["FP_FILTER2"]["TRANS"][:idl_filter_array_length]

    assert u.allclose(idl_filter_transmission_auto, filter_transmission)


@pytest.mark.parametrize("channel_name", channel_names)
def test_filter2_mesh_transmission(channel_name):
    channel_filter = Channel(channel_name)
    filter_mesh_transmission = channel_filter.filter_2.mesh_trans

    idl_filter_mesh_transmission_auto = v6_genx_s[
        _channel_name_to_index_mapping[channel_name]
    ]["FP_FILTER2"]["MESH_TRANS"]

    assert filter_mesh_transmission == idl_filter_mesh_transmission_auto


@pytest.mark.parametrize("channel_name", channel_names)
def test_filter2_substrate(channel_name):
    channel_filter = Channel(channel_name)
    filter_substrate = channel_filter.filter_2.substrate

    idl_filter_substrate_auto = v6_genx_s[_channel_name_to_index_mapping[channel_name]][
        "FP_FILTER2"
    ]["SUBSTRATE"]

    assert filter_substrate == idl_filter_substrate_auto


@pytest.mark.parametrize("channel_name", channel_names)
def test_geometry_name(channel_name):
    channel_filter = Channel(channel_name)
    geometry_name = channel_filter.geometry.name

    IDL_geometry_name_AUTO = v6_genx_s[_channel_name_to_index_mapping[channel_name]][
        "GEOM"
    ]["LONG_NAME"]

    assert geometry_name == IDL_geometry_name_AUTO


@pytest.mark.parametrize("channel_name", channel_names)
def test_geometry_focal_len(channel_name):
    channel_filter = Channel(channel_name)
    geometry_focal_len = channel_filter.geometry.focal_len

    IDL_geometry_focal_len_AUTO = (
        v6_genx_s[_channel_name_to_index_mapping[channel_name]]["GEOM"]["FOC_LEN"]
        * u.cm
    )

    assert u.isclose(geometry_focal_len, IDL_geometry_focal_len_AUTO)


@pytest.mark.parametrize("channel_name", channel_names)
def test_geometry_aperture_area(channel_name):
    channel_filter = Channel(channel_name)
    geometry_aperture_area = channel_filter.geometry.aperture_area

    idl_geometry_aperture_area_auto = (
        v6_genx_s[_channel_name_to_index_mapping[channel_name]]["GEOM"]["APERTURE_AREA"]
        * u.cm ** 2
    )

    assert u.isclose(geometry_aperture_area, idl_geometry_aperture_area_auto)


@pytest.mark.parametrize("channel_name", channel_names)
def test_mirror1_name(channel_name):
    channel_filter = Channel(channel_name)
    mirror_name = channel_filter.mirror_1.name

    IDL_mirror_name_AUTO = v6_genx_s[_channel_name_to_index_mapping[channel_name]][
        "MIRROR1"
    ]["LONG_NAME"]

    assert mirror_name == IDL_mirror_name_AUTO


@pytest.mark.parametrize("channel_name", channel_names)
def test_mirror1_material(channel_name):
    channel_filter = Channel(channel_name)
    mirror_material = channel_filter.mirror_1.material

    IDL_mirror_material_AUTO = v6_genx_s[_channel_name_to_index_mapping[channel_name]][
        "MIRROR1"
    ]["MATERIAL"]

    assert mirror_material == IDL_mirror_material_AUTO


@pytest.mark.parametrize("channel_name", channel_names)
def test_mirror1_density(channel_name):
    channel_filter = Channel(channel_name)
    mirror_density = channel_filter.mirror_1.density

    idl_mirror_density_auto = v6_genx_s[_channel_name_to_index_mapping[channel_name]][
        "MIRROR1"
    ]["DENS"] * (u.g * u.cm ** -3)

    assert u.isclose(mirror_density, idl_mirror_density_auto)


@pytest.mark.parametrize("channel_name", channel_names)
def test_mirro1_graze_angle(channel_name):
    channel_filter = Channel(channel_name)
    mirror_graze_angle = channel_filter.mirror_1.graze_angle

    idl_mirror_graze_angle_auto = (
        v6_genx_s[_channel_name_to_index_mapping[channel_name]]["MIRROR1"][
            "GRAZE_ANGLE"
        ]
        * u.deg
    )

    assert u.isclose(mirror_graze_angle, idl_mirror_graze_angle_auto)

    idl_mirror_graze_angle_manu = [0.910000] * u.deg
    assert u.isclose(idl_mirror_graze_angle_manu, mirror_graze_angle)


@pytest.mark.parametrize("channel_name", channel_names)
def test_mirror1_wavelength(channel_name):
    channel_filter = Channel(channel_name)

    mirror_number_of_length = int(channel_filter.mirror_1.number_of_wavelengths)
    mirror_wavelength = channel_filter.mirror_1.wavelength[:mirror_number_of_length]

    idl_mirror_array_length = int(
        v6_genx_s[_channel_name_to_index_mapping[channel_name]]["MIRROR1"]["LENGTH"]
    )
    idl_mirror_wavelength_auto = v6_genx_s[
        _channel_name_to_index_mapping[channel_name]
    ]["MIRROR1"]["WAVE"][:idl_mirror_array_length] * u.Unit("Angstrom")

    assert u.allclose(idl_mirror_wavelength_auto, mirror_wavelength)

    idl_mirror_wavelength_manu = [
        1.00000,
        1.10000,
        1.20000,
        1.30000,
        1.40000,
        1.50000,
        1.60000,
        1.70000,
        1.80000,
        1.90000,
    ] * u.angstrom
    assert u.allclose(idl_mirror_wavelength_manu, mirror_wavelength[0:10])


@pytest.mark.parametrize("channel_name", channel_names)
def test_mirror1_reflection(channel_name):
    channel_filter = Channel(channel_name)

    mirror_number_of_length = int(channel_filter.mirror_1.number_of_wavelengths)
    mirror_reflection = channel_filter.mirror_1.reflection[:mirror_number_of_length]

    idl_mirror_array_length = int(
        v6_genx_s[_channel_name_to_index_mapping[channel_name]]["MIRROR1"]["LENGTH"]
    )
    idl_mirror_wavelength_auto = (
        v6_genx_s[_channel_name_to_index_mapping[channel_name]]["MIRROR1"]["REFL"][
            :idl_mirror_array_length
        ]
        * u.angstrom
    )

    assert u.allclose(idl_mirror_wavelength_auto, mirror_reflection)


@pytest.mark.parametrize("channel_name", channel_names)
def test_mirror2_name(channel_name):
    channel_filter = Channel(channel_name)
    mirror_name = channel_filter.mirror_1.name

    idl_mirror_name_auto = v6_genx_s[_channel_name_to_index_mapping[channel_name]][
        "MIRROR2"
    ]["LONG_NAME"]

    assert mirror_name == idl_mirror_name_auto


@pytest.mark.parametrize("channel_name", channel_names)
def test_mirror2_material(channel_name):
    channel_filter = Channel(channel_name)
    mirror_material = channel_filter.mirror_1.material

    idl_mirror_material_auto = v6_genx_s[_channel_name_to_index_mapping[channel_name]][
        "MIRROR2"
    ]["MATERIAL"]

    assert mirror_material == idl_mirror_material_auto


@pytest.mark.parametrize("channel_name", channel_names)
def test_mirror2_density(channel_name):
    channel_filter = Channel(channel_name)
    mirror_density = channel_filter.mirror_1.density

    idl_mirror_density_auto = v6_genx_s[_channel_name_to_index_mapping[channel_name]][
        "MIRROR2"
    ]["DENS"] * (u.g * u.cm ** -3)

    assert u.isclose(mirror_density, idl_mirror_density_auto)


@pytest.mark.parametrize("channel_name", channel_names)
def test_mirror2_graze_angle(channel_name):
    channel_filter = Channel(channel_name)
    mirror_graze_angle = channel_filter.mirror_1.graze_angle

    idl_mirror_graze_angle_auto = (
        v6_genx_s[_channel_name_to_index_mapping[channel_name]]["MIRROR2"][
            "GRAZE_ANGLE"
        ]
        * u.deg
    )

    assert u.isclose(mirror_graze_angle, idl_mirror_graze_angle_auto)

    idl_mirror_graze_angle_manu = [0.910000] * u.deg
    assert u.isclose(idl_mirror_graze_angle_manu, mirror_graze_angle)


@pytest.mark.parametrize("channel_name", channel_names)
def test_mirror2_wavelength(channel_name):
    channel_filter = Channel(channel_name)

    mirror_number_of_length = int(channel_filter.mirror_1.number_of_wavelengths)
    mirror_wavelength = channel_filter.mirror_1.wavelength[:mirror_number_of_length]

    idl_mirror_array_length = int(
        v6_genx_s[_channel_name_to_index_mapping[channel_name]]["MIRROR2"]["LENGTH"]
    )
    idl_mirror_wavelength_auto = (
        v6_genx_s[_channel_name_to_index_mapping[channel_name]]["MIRROR2"]["WAVE"][
            :idl_mirror_array_length
        ]
        * u.angstrom
    )

    assert u.allclose(idl_mirror_wavelength_auto, mirror_wavelength)

    idl_mirror_wavelength_manu = [
        1.00000,
        1.10000,
        1.20000,
        1.30000,
        1.40000,
        1.50000,
        1.60000,
        1.70000,
        1.80000,
        1.90000,
    ] * u.angstrom
    assert u.allclose(idl_mirror_wavelength_manu, mirror_wavelength[0:10])


@pytest.mark.parametrize("channel_name", channel_names)
def test_mirror2_reflection(channel_name):
    channel_filter = Channel(channel_name)

    mirror_number_of_length = int(channel_filter.mirror_1.number_of_wavelengths)
    mirror_reflection = channel_filter.mirror_1.reflection[:mirror_number_of_length]

    idl_mirror_array_length = int(
        v6_genx_s[_channel_name_to_index_mapping[channel_name]]["MIRROR2"]["LENGTH"]
    )
    idl_mirror_wavelength_auto = (
        v6_genx_s[_channel_name_to_index_mapping[channel_name]]["MIRROR2"]["REFL"][
            :idl_mirror_array_length
        ]
        * u.angstrom
    )

    assert u.allclose(idl_mirror_wavelength_auto, mirror_reflection)


@pytest.mark.parametrize("channel_name", channel_names)
def test_channel_name(channel_name):
    channel_filter = Channel(channel_name)
    name = channel_filter.name

    IDL_mirror_name_AUTO = v6_genx_s[_channel_name_to_index_mapping[channel_name]][
        "NAME"
    ]

    assert name == IDL_mirror_name_AUTO


@pytest.mark.parametrize("channel_name", channel_names)
def test_channel_wavelength(channel_name):
    channel_filter = Channel(channel_name)

    wavelength_length = int(channel_filter.number_of_wavelengths)
    wavelength = channel_filter.wavelength[:wavelength_length]

    idl_array_length = int(
        v6_genx_s[_channel_name_to_index_mapping[channel_name]]["LENGTH"]
    )
    idl_wavelength_auto = (
        v6_genx_s[_channel_name_to_index_mapping[channel_name]]["WAVE"][
            :idl_array_length
        ]
        * u.angstrom
    )

    assert u.allclose(idl_wavelength_auto, wavelength)

    idl_mirror_wavelength_manu = [
        9.00000,
        9.10000,
        9.20000,
        9.30000,
        9.40000,
        9.50000,
        9.60000,
        9.70000,
        9.80000,
        9.90000,
    ] * u.angstrom
    assert u.allclose(idl_mirror_wavelength_manu, wavelength[80:90])


@pytest.mark.parametrize("channel_name", channel_names)
def test_channel_transmission(channel_name):
    channel_filter = Channel(channel_name)

    transmission_length = int(channel_filter.number_of_wavelengths)
    transmission = channel_filter.transmission[:transmission_length]

    idl_array_length = int(
        v6_genx_s[_channel_name_to_index_mapping[channel_name]]["LENGTH"]
    )
    idl_transmission_auto = v6_genx_s[_channel_name_to_index_mapping[channel_name]][
        "TRANS"
    ][:idl_array_length]

    assert u.allclose(idl_transmission_auto, transmission)


@pytest.mark.parametrize("channel_name", channel_names)
def test_channel_number_of_wavelengths(channel_name):
    channel_filter = Channel(channel_name)

    channel_number_of_wavelengths = channel_filter.number_of_wavelengths

    idl_array_length = v6_genx_s[_channel_name_to_index_mapping[channel_name]]["LENGTH"]
    assert channel_number_of_wavelengths == idl_array_length


@pytest.mark.parametrize("channel_name", channel_names)
def test_channel_observatory(channel_name):
    channel_filter = Channel(channel_name)

    observatory = channel_filter.observatory

    idl_observatory = v6_genx_s[_channel_name_to_index_mapping[channel_name]][
        "OBSERVATORY"
    ]

    assert observatory == idl_observatory


@pytest.mark.parametrize("channel_name", channel_names)
def test_channel_instrument(channel_name):
    channel_filter = Channel(channel_name)

    instrument = channel_filter.instrument

    idl_instrument = v6_genx_s[_channel_name_to_index_mapping[channel_name]][
        "INSTRUMENT"
    ]

    assert instrument == idl_instrument


@pytest.mark.parametrize("attr", ["wavelength", "number_of_wavelengths"])
def test_open_channel(attr):

    open_filter = Channel("open")
    sample_filter = Channel("Al-mesh")

    open_value = getattr(open_filter, attr)
    sample_value = getattr(sample_filter, attr)
    assert u.allclose(open_value, sample_value)


def test_open_transmission():
    open_filter = Channel("open")
    assert u.allclose(open_filter.transmission, 1)
