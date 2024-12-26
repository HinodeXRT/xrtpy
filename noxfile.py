import nox

nox.options.sessions = ["tests"]
python_versions = ("3.10", "3.11", "3.12")


@nox.session
def tests(session):
    """
    Run tests with pytest.
    """
    pytest_options = {}
    session.install(".[tests]")
    session.run("pytest", *pytest_options)


@nox.session
def linters(session):
    """
    Run all pre-commit hooks on all files.
    """
    session.install("pre-commit")
    session.run("pre-commit", "run", "--all-files", *session.posargs)


@nox.session
def import_package(session):
    """
    Import xrtpy.
    """
    session.install(".")
    session.run("python", "-c", "import xrtpy")


@nox.session
def docs(session):
    """
    Build documentation with Sphinx.
    """
    sphinx_paths = ["docs", "docs/_build/html"]
    sphinx_fail_on_warnings = ["-W", "--keep-going"]
    sphinx_builder = ["-b", "html"]
    sphinx_nitpicky = ["-n"]
    sphinx_opts = (
        sphinx_paths + sphinx_fail_on_warnings + sphinx_builder + sphinx_nitpicky
    )
    session.install(".[docs]")
    session.run(
        "sphinx-build",
        *sphinx_opts,
        *sphinx_nitpicky,
        *session.posargs,
    )
