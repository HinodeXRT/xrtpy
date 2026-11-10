## Instructions

**If all checks pass ✅, please merge this pull request.** If any checks fail due to a breaking change in a dependency 🚨, please address the problems before merging.

## Description

This pull request upgrades [`uv.lock`]: the cross-platform lockfile that specifies the Python environments used when running tests, building documentation, and performing continuous integration (CI) checks. Locking and periodically updating the Python environment lets us quarantine breaking changes before they start causing spontaneous failures on unrelated pull requests, while ensuring that everyone is using the same versions of dependencies to perform checks.

The lockfile was upgraded via [`uv lock --upgrade`]. This workflow runs the [ty] static type checker to update `# ty:ignore` comments and perform autofixes, following by any autofixes made when running [pre-commit].

> [!NOTE]
> If it is necessary to temporarily place an upper limit on a dependency in [`pyproject.toml`], please [create an issue] to remove this upper limit before the next release.

> [!CAUTION]
> This workflow contains an experimental step to [bump-minimum-requirements].

## CI snapshot

[![CI](https://github.com/PlasmaPy/PlasmaPy/actions/workflows/ci.yml/badge.svg)](https://github.com/PlasmaPy/PlasmaPy/actions/workflows/ci.yml)
[![comprehensive tests](https://github.com/PlasmaPy/PlasmaPy/actions/workflows/ci-comprehensive.yml/badge.svg)](https://github.com/PlasmaPy/PlasmaPy/actions/workflows/ci-comprehensive.yml)
[![upstream tests](https://github.com/PlasmaPy/PlasmaPy/actions/workflows/upstream-tests.yml/badge.svg)](https://github.com/PlasmaPy/PlasmaPy/actions/workflows/upstream-tests.yml)
[![upstream docs](https://github.com/PlasmaPy/PlasmaPy/actions/workflows/upstream-docs.yml/badge.svg)](https://github.com/PlasmaPy/PlasmaPy/actions/workflows/upstream-docs.yml)
[![PyHC Actions](https://github.com/PlasmaPy/PlasmaPy/actions/workflows/pyhc-actions.yml/badge.svg)](https://github.com/PlasmaPy/PlasmaPy/actions/workflows/pyhc-actions.yml)
[![upgrade uv.lock](https://github.com/PlasmaPy/PlasmaPy/actions/workflows/upgrade-uv-lock.yml/badge.svg)](https://github.com/PlasmaPy/PlasmaPy/actions/workflows/upgrade-uv-lock.yml)

## Dependencies upgraded by workflow

[bump-minimum-requirements]: https://github.com/namurphy/bump-minimum-requirements
[create an issue]: https://github.com/PlasmaPy/PlasmaPy/issues/new?title=Remove+upper+limit+on+version+of
[pre-commit]: https://pre-commit.com
[ty]: https://docs.astral.sh/ty
[`pyproject.toml`]: https://github.com/PlasmaPy/PlasmaPy/blob/main/pyproject.toml
[`uv lock --upgrade`]: https://docs.astral.sh/uv/reference/cli/#uv-lock
[`uv.lock`]: https://docs.astral.sh/uv/guides/projects/#uvlock
