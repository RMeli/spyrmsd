# Development, testing, and deployment tools

This directory contains a collection of tools for running Continuous Integration (CI) tests, 
conda installation, and other development tools not directly related to the coding process.

## Continuous Integration

The code is tested on Ubuntu and macOS ([Travis-CI](https://about.travis-ci.com/)) and Windows ([AppVeyor](https://www.appveyor.com/)).

Content:
* `travis-ci`: Linux and macOS based testing through  
  * `before_install.sh`: Miniconda pre-package installation script for Travis

## Conda Environments

Conda is the recommended package manager for this project.

* `conda-envs`: directory containing the YAML file(s) which fully describe Conda environments
  * `spyrmsd.yaml`: Full Conda environment for `spyrmsd`.
  
Channels are usually not specified here and therefore respect global Conda configuration.

### Scripts

* `scripts`
  * `create_conda_env.py`: Helper program for spinning up new Conda environments based on a starter file with Python Version and Env.

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
