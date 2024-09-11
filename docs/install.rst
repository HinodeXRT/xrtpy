.. _xrtpy-install:

****************
Installing XRTpy
****************

Installing Python
=================

There are many ways to install Python, but even if you have Python installed somewhere on your computer we recommend following these instructions anyway.
That's because we will create a new Python environment.
As well as containing a Python installation, this environment provides an isolated place to install Python packages (like ``xrtpy``) without affecting any other current Python installation.
If you already have Python and ``conda`` working you can skip the next section.
Please note that XRTpy requires Python_ |minpython| or newer.

.. tip::

   New versions of Python_ are released annually in October, and it can take a few months for the scientific Python ecosystem to catch up.
   If you have trouble installing `xrtpy` on the most recent Python_ version between October and ∼March, then try installing it on the second most recent version.

`If you are using Anaconda, we recommend that you uninstall it as the default package channel(s) have a restrictive license which means you might not be able to use it for free <https://sunpy.org/posts/2024/2024-08-09-anaconda/>`__.
Instead, we recommend that you use miniforge which is a minimal installer that setups conda with the conda-forge channel, which is free to use for everyone.
If you are using miniforge, you can skip the next section

Installing miniforge
--------------------

If you don't already have a Python installation then we recommend installing Python with `miniforge <https://github.com/conda-forge/miniforge/#miniforge>`__.
This will install ``conda`` and automatically configure the default channel (a channel is a remote software repository) to be ``conda-forge``, which is where ``sunpy`` is available.

First, download the installer for your system and architecture from the links below:

.. grid:: 3

    .. grid-item-card:: Linux

        `x86-64 <https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh>`__

        `aarch64 <https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-aarch64.sh>`__

        `ppc64le <https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-ppc64le.sh>`__

    .. grid-item-card:: Windows
        :link: https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Windows-x86_64.exe

        `x86-64 <https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Windows-x86_64.exe>`__

    .. grid-item-card:: Mac

        `x86-64 <https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-x86_64.sh>`__

        `arm64 (Apple
        Silicon) <https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-arm64.sh>`__

Then select your platform to install miniforge:

.. tab-set::

    .. tab-item:: Linux & Mac
        :sync: platform

        Linux & Mac Run the script downloaded above, with
        ``bash <filename>``. The following should work:

        .. code-block:: console

            bash Miniforge3-$(uname)-$(uname -m).sh

        Once the installer has completed, close and reopen your terminal.

    .. tab-item:: Windows
        :sync: platform

        Double click the executable file downloaded from
        the links above.

        Once the installer has completed you should have a new "miniforge
        Prompt" entry in your start menu.

In a new terminal (miniforge Prompt on Windows) run ``conda list`` to test that the install has worked.

Installing xrtpy
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

Once the environment active, to acquire a full ``xrtpy`` installation:

.. code-block:: bash

    $ pip install xrtpy


.. warning::

    If you get a ``PermissionError`` this means that you do not have the required administrative access to install new packages to your Python installation.
    Do **not** install any Python packages using ``sudo``.
    This error implies you have an incorrectly configured virtual environment or it is not activated.

.. note::

   If you noticed any places where the installation instructions could be improved or have become out of date, please create an issue on `XRTpy's GitHub repository`_.
   It would really help!
