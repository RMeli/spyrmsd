# Development, testing, and deployment tools

This directory contains a collection of tools for running Continuous Integration (CI) tests,
conda installation, and other development tools not directly related to the coding process.

## Continuous Integration

The code is tested on Ubuntu, macOS and Windows ([GitHub Actions](https://docs.github.com/en/actions)).

## Conda Environments

Conda is the recommended package manager for this project.

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
