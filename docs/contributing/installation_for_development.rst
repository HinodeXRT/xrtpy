.. _installation for development:

****************************
Installation for Development
****************************

Development Environment
=======================
To set up your development environment:

1. Clone the repository::
   .. code-block:: shell
      git clone https://github.com/HinodeXRT/xrtpy.git
      cd xrtpy

2. Install the package and required dependencies::
   .. code-block:: shell
      pip install -e .[dev,docs,tests]

Coding Standards
================
- Follow the `PEP-8`_ coding style.
- Write clear and concise commit messages.
- Include docstrings for all functions and classes.
- Ensure that your code is covered by tests and that all tests pass before submitting a PR.

Testing
=======
We use `pytest` for testing, with Nox_ as the test runner. To run the
tests locally, use the following command in the top-level directory:

.. code-block:: shell
   nox

Ensure that all tests pass before merging your PR.

Documentation
=============
We use Sphinx_ to build documentation via a Nox_ session. To build
documentation locally, run

.. code-block:: shell
   nox -s docs

Communication
=============
For any questions or discussions, you can email us at `xrtpy@cfa.harvard.edu`.

.. _PEP-8: https://peps.python.org/pep-0008
.. _Nox: https://nox.thea.codes
.. _Sphinx: https://www.sphinx-doc.org