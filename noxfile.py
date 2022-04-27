import nox

python_versions = ("3.8", "3.9", "3.10")


@nox.session(python=python_versions)
def tests(session):
    session.install("-r", "requirements/tests.txt")
    session.install(".")
    session.run("pytest")


@nox.session
def linters(session):
    session.install("-r", "requirements/tests.txt")
    flake8_options = ["--count", "--show-source", "--statistics"]
    session.run("flake8", "plasmapy", *flake8_options, *session.posargs)


@nox.session
def build_docs(session):
    session.install("-r", "requirements/docs.txt")
    session.install(".")
    session.run(
        "sphinx-build",
        "docs",
        "docs/_build/html",
        "-W",
        "--keep-going",
        "-b",
        "html",
        *session.posargs,
    )


@nox.session
def build_docs_no_examples(session):
    session.install("-r", "requirements/docs.txt")
    session.install(".")
    session.run(
        "sphinx-build",
        "-D",
        "nbsphinx_execute=never",
        "docs",
        "docs/_build/html",
        "-W",
        "--keep-going",
        "-b",
        "html",
        *session.posargs,
    )


@nox.session(python="3.10")
def codespell(session):
    session.install("codespell")
    session.run("codespell", ".")
