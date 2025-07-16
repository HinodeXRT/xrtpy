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


def validate_and_format_filters(filters):
    """
    Validates and formats XRT filter names. Ensures no duplicates.

    Parameters
    ----------
    filters : str or list of str
        One or more filter names provided by the user.

    Returns
    -------
    list of str
        A list of unique, standardized filter names.

    Raises
    ------
    TypeError
        If input is not a string or list of strings.
    ValueError
        If duplicate filters are found after formatting.

    Examples
    --------
    >>> validate_and_format_filters("al_poly")
    ['Al-Poly']

    >>> validate_and_format_filters(["be_thin", "be_thin"])
    ValueError: Duplicate filters detected: ['Be-Thin', 'Be-Thin']
    """
    if isinstance(filters, str):
        filters = [filters]
    elif not isinstance(filters, list):
        raise TypeError("Input must be a string or list of strings.")
    
    if not all(isinstance(f, str) for f in filters):
        raise TypeError("All filter names must be strings.")

    formatted = [solve_filter_name(f) for f in filters]

    # Check for duplicates
    if len(set(formatted)) != len(formatted):
        raise ValueError(f"Duplicate filters detected: {formatted}")

    return formatted