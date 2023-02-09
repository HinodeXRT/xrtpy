.. |bibliography| replace:: :ref:`bibliography`\
.. |glossary| replace:: :ref:`glossary`\
.. |inf| replace:: `~numpy.inf`
.. |Map| replace:: `sunpy.map.Map`
.. |minpython| replace:: 3.8
.. |nan| replace:: `~numpy.nan`
.. |ndarray| replace:: :class:`~numpy.ndarray`
.. |Quantity| replace:: :class:`~astropy.units.Quantity`
.. |TimeDelta| replace:: :class:`~astropy.time.TimeDelta`
.. |Time| replace:: :class:`~astropy.time.Time`
.. |Channel| replace:: :class:`~xrtpy.response.channel.Channel`
.. |EffectiveAreaFundamental| replace:: :class:`~xrtpy.response.effective_area.EffectiveAreaFundamental`
.. |TemperatureResponseFundamental| replace:: :class:`~xrtpy.response.temperature_response.TemperatureResponseFundamental`

.. A workaround for nested inline literals so that the filename will get
   formatted like a file but will be a link. In the text, these get used
   with the syntax for a substitution followed by an underscore to
   indicate that it's for a link: |docs/_static|_

.. For these workarounds, if the replacement is something in single back
   ticks (e.g., `xarray`), then it should also be added to
   nitpick_ignore_regex in docs/conf.py so that it doesn't get counted
   as an error in a nitpicky doc build (e.g., tox -e doc_build_nitpicky).

.. _`docs/_static`: https://github.com/HinodeXRT/xrtpy/tree/main/docs/_static
.. |docs/_static| replace:: :file:`docs/_static`

.. _`docs/api_static`: https://github.com/HinodeXRT/xrtpy/tree/main/docs/api_static
.. |docs/api_static| replace:: :file:`docs/api_static`

.. _`docs/conf.py`: https://github.com/HinodeXRT/xrtpy/blob/main/docs/conf.py
.. |docs/conf.py| replace:: :file:`docs/conf.py`

.. _`docs/glossary.rst`: https://github.com/HinodeXRT/PlasmaPy/blob/main/docs/glossary.rst
.. |docs/glossary.rst| replace:: :file:`docs/glossary.rst`

.. _`docs/bibliography.bib`: https://github.com/HinodeXRT/xrtpy/blob/main/docs/bibliography.bib
.. |docs/bibliography.bib| replace:: :file:`docs/bibliography.bib`

.. _`setup.cfg`: https://github.com/HinodeXRT/xrtpy/blob/main/setup.cfg
.. |setup.cfg| replace:: :file:`setup.cfg`
