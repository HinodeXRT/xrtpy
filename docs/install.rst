.. _XRTpy-install:

****************
Installing XRTpy
****************

Installing Python
=================

To use ``xrtpy``, you will need to have Python installed on your system.
We recommend following :ref:`sunpy-tutorial-installing` which walks you through installing Python.

Installing XRTpy
----------------

To install ``xrtpy``, start by launching a terminal (under a UNIX-like system) or the miniforge Prompt (under Windows).
Now we will create and activate a new virtual environment to install ``xrtpy`` into:

.. code-block:: bash

    $ conda create --name xrtpy
    # Only run the following two lines
    # if you have NOT installed miniforge or added conda-forge to your channels
    # Do not run these lines if you are using Anaconda
    $ conda config --add channels conda-forge
    $ conda config --set channel_priority strict
    $ conda activate xrtpy

In this case the environment is named 'xrtpy'.
Feel free to change this to a different environment name.

The benefit of using a virtual environment is that it allows you to install packages without affecting any other Python installation on your system.
This also means you can work on multiple projects (research or coding) with different package requirements without them interfering with each other.

Now we have a fresh environment we can install ``xrtpy``:

.. code-block:: bash

    $ conda install xrtpy

This will install ``xrtpy`` and all of its dependencies.
If you want to install another package later, you can run ``conda install <package_name>``.

pip
~~~

This is for installing ``xrtpy`` within a Python environment, where ``pip`` has been used to install all previous packages.
You will want to make sure you are using a `Python virtual environment <https://packaging.python.org/guides/installing-using-pip-and-virtual-environments/>`__.

Once the environment active, to acquire a ``xrtpy`` installation:

.. code-block:: bash

    $ pip install xrtpy


.. warning::

    If you get a ``PermissionError`` this means that you do not have the required administrative access to install new packages to your Python installation.
    Do **not** install any Python packages using ``sudo``.
    This error implies you have an incorrectly configured virtual environment or it is not activated.

.. note::

   If you noticed any places where the installation instructions could be improved or have become out of date, please create an issue on `XRTpy's GitHub repository`_.
   It would really help!
