import sys

import nox

nox.options.default_venv_backend = "uv"

supported_python_versions = ("3.11", "3.12", "3.13")

maxpython = max(supported_python_versions)
minpython = min(supported_python_versions)
docpython = "3.13"

current_python = f"{sys.version_info.major}.{sys.version_info.minor}"

nox.options.sessions = [f"tests-{current_python}(all)"]

pytest_command: tuple[str, ...] = (
    "pytest",
    "--pyargs",
    "--durations=5",
    "--tb=short",
    "-n=auto",
    "--dist=loadfile",
)

with_coverage: tuple[str, ...] = (
    "--cov=plasmapy",
    "--cov-report=xml",
    "--cov-config=pyproject.toml",
    "--cov-append",
    "--cov-report",
    "xml:coverage.xml",
)

test_specifiers: list = [
    nox.param("run all tests", id="all"),
    nox.param("with code coverage", id="cov"),
    nox.param("lowest-direct", id="lowest-direct"),
]


@nox.session(python=supported_python_versions)
@nox.parametrize("test_specifier", test_specifiers)
def tests(session, test_specifier: nox._parametrize.Param) -> None:
    """Run tests with pytest."""

    install_options = (
        ["--resolution=lowest-direct"] if test_specifier == "lowest-direct" else []
    )

    pytest_options: list[str] = (
        with_coverage if test_specifier == "with code coverage" else []
    )

    session.install("uv")
    session.install(".[tests]", *install_options)

    session.run("pytest", *pytest_options, *session.posargs)


@nox.session
def lint(session: nox.Session) -> None:
    """Run all pre-commit hooks on all files."""
    session.install("pre-commit")
    session.run(
        "pre-commit",
        "run",
        "--all-files",
        "--show-diff-on-failure",
        *session.posargs,
    )


@nox.session
def build(session: nox.Session) -> None:
    """Build & verify the source distribution and wheel."""
    session.install("twine", "build")
    build_command = ("python", "-m", "build")
    session.run(*build_command, "--sdist")
    session.run(*build_command, "--wheel")
    session.run("twine", "check", "dist/*", *session.posargs)


@nox.session(python=docpython)
def docs(session):
    """
    Build documentation with Sphinx.

    This session may require installation of pandoc and graphviz.
    """
    sphinx_paths = ["docs", "docs/_build/html"]
    sphinx_fail_on_warnings = ["-W", "--keep-going"]
    sphinx_builder = ["-b", "html"]
    sphinx_nitpicky = ["-n"]
    sphinx_opts = (
        sphinx_paths + sphinx_fail_on_warnings + sphinx_builder + sphinx_nitpicky
    )
    session.install(".[docs]", "--exclude-newer=2025-03-01")
    session.run(
        "sphinx-build",
        *sphinx_opts,
        *sphinx_nitpicky,
        *session.posargs,
    )
