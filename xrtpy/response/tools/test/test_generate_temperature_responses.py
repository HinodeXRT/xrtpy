import astropy.units as u
import pytest
from astropy.time import Time

from xrtpy.response.tools import generate_temperature_responses


@pytest.mark.parametrize(
    "filters",
    [
        (["Al-mesh"]),
        (["Al-poly"]),
        (["Al-mesh", "Be-thick"]),
        (["Al-mesh", "Be-thin", "C-poly", "Ti-poly"]),
    ],
)
def test_response_length_matches_filters(filters):
    responses = generate_temperature_responses(
        filters, "2011-01-28T11:02:31", "Photospheric"
    )
    assert len(responses) == len(filters)


@pytest.mark.parametrize(
    "obs_date",
    [
        "2008-01-01T00:00:00",
        "2011-01-28T11:02:31",
        "2022-06-15T12:30:00",
    ],
)
def test_responses_for_different_dates(obs_date):
    filters = ["Al-mesh", "Be-thin", "C-poly", "Ti-poly"]
    responses = generate_temperature_responses(filters, obs_date, "Photospheric")
    assert len(responses) == len(filters)
    for r in responses:
        assert isinstance(r.temperature, u.Quantity)
        assert isinstance(r.response, u.Quantity)


@pytest.mark.parametrize("abundance", ["Coronal", "Photospheric", "Hybrid"])
def test_responses_for_different_abundance_models(abundance):
    filters = ["Al-mesh", "Be-thick"]
    responses = generate_temperature_responses(
        filters, "2011-01-28T11:02:31", abundance
    )
    assert len(responses) == len(filters)
    for r in responses:
        assert r.temperature.unit.is_equivalent(u.K)
        assert r.response.unit.is_equivalent(u.DN * u.cm**5 / (u.pix * u.s))


def test_shape_and_units_consistency():
    filters = ["Al-poly", "Ti-poly"]
    responses = generate_temperature_responses(
        filters, "2011-01-28T11:02:31", "Coronal"
    )
    for r in responses:
        assert len(r.temperature) == len(r.response)
        assert r.temperature.unit.is_equivalent(u.K)
        assert r.response.unit.is_equivalent(u.DN * u.cm**5 / (u.pix * u.s))


def test_astropy_time_input():
    filters = ["Al-mesh"]
    time_obj = Time("2011-01-28T11:02:31")
    responses = generate_temperature_responses(filters, time_obj, "Coronal")
    assert len(responses) == len(filters)


# def test_invalid_filter_raises_exception():
#     """
#     Test that passing an invalid filter name raises a ValueError.
#     """
#     filters = ["Not-a-real-filter"]
#     with pytest.raises(ValueError, match="Invalid filter name"):
#         generate_temperature_responses(filters, "2011-01-28T11:02:31", "Coronal")


def test_invalid_filter_raises_exception():
    """
    Test that passing an invalid filter name raises an exception.
    """
    filters = ["Not-a-real-filter"]
    with pytest.raises(Exception):
        generate_temperature_responses(filters, "2011-01-28T11:02:31", "Coronal")
