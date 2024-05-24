.. _contributing:

*********************
Contributing to XRTpy
*********************

Thank you for your interest in contributing to XRTpy!
We welcome contributions from the community to improve and expand the functionality of this package.

There are several ways you can contribute to XRTpy:

1. **Reporting Issues**: If you encounter any bugs or have suggestions for improvements, please report them using the `GitHub-issue`_ page. Provide as much detail as possible, including steps to reproduce the issue and any relevant screenshots or code snippets.

2. **Submitting Pull Requests**: If you want to contribute code, follow these steps:
   - Fork the repository on GitHub.
   - Create a new branch from the `main` branch for your changes.
   - Make your changes, ensuring that you follow the coding standards and include tests for any new functionality.
   - Commit your changes with clear and descriptive commit messages.
   - Push your branch to your forked repository.
   - Create a pull request (PR) from your branch to the `main` branch of the original repository.

3. **Improving Documentation**: Good documentation is crucial for the usability of the package. You can help by improving existing documentation, or writing new tutorials.


Development Environment
========================
To set up your development environment:

1. Clone the repository::

   ```git clone https://github.com/HinodeXRT/xrtpy.git```
   `cd xrtpy`

2. Install the required dependencies::

   `pip install -r requirements.txt1

3. Install the package in editable mode::

   `pip install -e .`

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


Thank you for contributing to XRTpy!


.. _PEP-8: https://peps.python.org/pep-0008/
.. _GitHub-issue: https://github.com/HinodeXRT/xrtpy/issues
