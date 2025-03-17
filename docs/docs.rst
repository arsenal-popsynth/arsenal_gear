============================
Arsenal Sphinx Documentation
============================
Arsenal gear uses the `Sphinx <http://sphinx-doc.org/>`_ documentation engine
to generate its documentation. The documentation is written in `reStructuredText
<http://docutils.sourceforge.net/rst.html>`_ (reST) and is located in the
``docs/`` directory of the `repository <https://github.com/arsenal-popsynth/arsenal_gear>`_. 

Testing the Documentation Locally
=================================

Before deploying a change to the documentation, you may want to test the build
locally. This is a fairly straightforward, four-step process.

1. Make sure you've actually installed arsenal gear.  From the root directory of the repository, run ``pip install .[all]``
2. Use the Makefile in the ``docs/`` directory to build the documentation: ``cd docs/ && make html`` 
3. Spin up an HTTP server in the directory containing the compiled documentation: ``cd _build/html && python -m http.server 1337`` 
4. Open a web browser and navigate to ``http://localhost:1337`` to view the documentation.

Deploying Documentation to GitHub Pages
=======================================

Arsenal gear uses GitHub actions to automatically deploy the documentation to `arsenal-popsynth.github.io/arsenal_gear <https://arsenal-popsynth.github.io/arsenal_gear>`_ whenever a change is merged to the ``main`` branch.
The actions workflow is defined in ``.github/workflows/docs.yml``. 

The workflow is triggered by a push to the ``main`` branch, and will ultimately put the html documents produced by Sphinx in the ``gh-pages`` branch of the repository. 

GitHub Pages is configured to serve the documentation from the ``gh-pages`` branch automatically, so you shouldn't need to do anything manually to deploy the documentation.