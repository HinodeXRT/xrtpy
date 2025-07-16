"""
Utility functions for validating and standardizing Hinode/XRT filter names.
"""

def solve_filter_name(name):
    """
    Standardizes an XRT filter name to match expected format.

    Parameters
    ----------
    name : str
        The filter name provided by the user. Can include hyphens, underscores, or inconsistent casing.

    Returns
    -------
    str
        A standardized filter name with correct capitalization and formatting.

    Raises
    ------
    TypeError
        If the provided name is not a string.

    Examples
    --------
        solve_filter_name("al_poly")
    'Al-Poly'

        solve_filter_name("be_thin/al_mesh")
    'Be-Thin/Al-Mesh'
    """
    if not isinstance(name, str):
        raise TypeError("name must be a string")

    name = name.replace("_", "-")
    parts = name.split("/")
    new_parts = [part.capitalize() for part in parts]
    return "/".join(new_parts)
