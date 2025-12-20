# Development, testing, and deployment tools

This directory contains a collection of tools for running Continuous Integration (CI) tests,
conda installation, and other development tools not directly related to the coding process.

## Continuous Integration

The code is tested on Ubuntu, macOS and Windows ([GitHub Actions](https://docs.github.com/en/actions)).

CI primarily uses [uv](https://github.com/astral-sh/uv) for fast and reliable dependency management. For configurations requiring graph-tool (which is only available via conda), micromamba is used in a hybrid approach.

## Conda Environments

Conda environments are provided for local development as an alternative to pip/uv.

* `conda-envs`: directory containing the YAML file(s) which fully describe Conda environments
  * `spyrmsd.yaml`: Full Conda environment for `spyrmsd`.

Channels are usually not specified here and therefore respect global Conda configuration.

## Deployment

### PyPI

Build wheel and sdist:

```bash
flit build
```

Upload wheel and sdist on [PyPI](https://pypi.org/):

```bash
flit publish
```
