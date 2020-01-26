# Compiling spyrmsd's Documentation

The docs for this project are built with [Sphinx](http://www.sphinx-doc.org/en/master/).

## Installation

Ensure that Sphinx and the ReadTheDocs theme are installed:

```bash
conda install sphinx sphinx_rtd_theme
```

## Documentation

### Automatic API Documentation

Create API documentation automatically with `sphinx-apidoc`:

```bash
sphinx-apidoc -f -M -e  -T -o source/api ../spyrmsd
```

### Build Documentation

Use the `Makefile` to compile static HTML pages:

```bash
make html
```

### Visualize Documentation

The compiled docs will be in the `build` directory and can be viewed by opening `index.html` (which may itself be inside a directory called `html/` depending on what version of Sphinx is installed).