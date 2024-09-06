.. _installation for development:

****************************
Installation for Development
****************************

Development Environment
========================
To set up your development environment:

1. Clone the repository::

      git clone https://github.com/HinodeXRT/xrtpy.git`
      cd xrtpy

2. Install the package in editable mode::

      pip install -e .

Coding Standards
================
- Follow the `PEP-8`_ coding style.
- Write clear and concise commit messages.
- Include docstrings for all functions and classes.
- Ensure that your code is covered by tests and that all tests pass before submitting a PR.

Testing
=======
We use `pytest` for testing. To run the tests, use the following command::

   pytest

Ensure that all tests pass before submitting your PR.

Communication
=============
For any questions or discussions, you can email us at `xrtpy@cfa.harvard.edu`.

.. _PEP-8: https://peps.python.org/pep-0008/
