History, status, and future
===========================
The `original PFSS implementation in Python <https://doi.org/10.5281/zenodo.1472183>`__ was written by `Anthony Yeates <https://www.maths.dur.ac.uk/users/anthony.yeates/>`__.
I (David Stansby) then took this and added tests, documentation, and integrated it with the SunPy project ecosystem to create the pfsspy package.
Most of the package was written by myself, with some community contributions.

I am the sole maintainer of the package, which is currently feature complete (providing a spherical PFSS solver in Python that integrates well with the SunPy package ecosystem).

Going forward the only changes to the package will be:
- Fixing `identified bugs<https://github.com/dstansby/pfsspy/issues?q=is%3Aopen+is%3Aissue+label%3ABug>`__
- Adding tightly scoped `new features<https://github.com/dstansby/pfsspy/issues?q=is%3Aopen+is%3Aissue+label%3AEnhancement>`__ to improve usability or integrate better with other SunPy packages.
- Improving the documentation.
- Retaining compatibility with future versions of the package dependencies (e.g. sunpy, Matplotlib...)
I intend to remain the sole maintainer of the package, but since I'm not currently active in solar physics research don't have lots of time to work on pfsspy. As a priority I will maintain the issue list with well scoped peices of work that others can contribute to fixing.

Community contributions
-----------------------
The `pfsspy issue tracker <https://github.com/dstansby/pfsspy/issues>`__ maintains a list of items that need working on.
Please open a new issue if you find a problem or want to see a new feature added!
If you want to fix any of the issues yourself, please open a pull request.
I will try my best to review and merge pull requests, but **please keep them as small as possible** to make my life easier!
General guidelines for contributing to open source can be found in the `sunpy developers guide <https://docs.sunpy.org/en/latest/dev_guide/index.html>`__
