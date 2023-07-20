What is climdex.pcic?
=====================

* `climdex.pcic` is a library for computing the `27 core indices of extreme climate`_. It was written for the `R statistical programming language`_ by the `Pacific Climate Impacts Consortium`_.

.. _27 core indices of extreme climate: http://etccdi.pacificclimate.org/list_27_indices.shtml
.. _R statistical programming language: http://www.r-project.org/
.. _Pacific Climate Impacts Consortium: http://pacificclimate.org/

Installation
============

`climdex.pcic` is not available from CRAN at this time.

One can install from the source tarball of our official release from
PCIC's website:

    # ensure dependency is installed
    > install.packages('PCICt')
    > install.packages('https://pacificclimate.org/R/climdex.pcic_1.1-11.tar.gz')

Or using the `devtools` package, one can install any arbitrary
version, commit or branch from GitHub:

    # if devtools is not installed
    > install.packages('devtools')

    # Install the master branch for bleeding edge code
    > devtools::install_github("pacificclimate/climdex.pcic")


Getting Help
============

New to programming or to R?
---------------------------

* Read the the `Software Carpentry`_  `Programming with R`_ lessons
* Read one of the man `R Manuals`_.
* Attend an `R Users Group`_ meeting.

.. _Software Carpentry: http://software-carpentry.org/index.html
.. _Programming with R: http://swcarpentry.github.io/r-novice-inflammation/
.. _R Manuals: http://cran.r-project.org/manuals.html
.. _R Users Group: http://r-users-group.meetup.com/

Looking for code?
-----------------

* Explore the `development repository`_ on GitHub.

.. _development repository: https://github.com/pacificclimate/climdex.pcic/

Need help using the package?
----------------------------

* Read the manual ::

    > library(climdex.pcic)
    Loading required package: PCICt
    > ?climdex.pcic

* Create a `new issue`_ on the `package issue tracker`_ and label it "help wanted"[1]_.

.. _new issue: https://github.com/pacificclimate/climdex.pcic/issues/new

Want to contribute?
-------------------

* To report a bug in pcic.climdex use the `package issue tracker`_ (after you've read the `bug reporting guide`_).
* To help with development read through the `contributor's guide`_

.. _bug reporting guide: https://github.com/pacificclimate/climdex.pcic/blob/master/CONTRIBUTING.rst#bug-reports
.. _package issue tracker: https://github.com/pacificclimate/climdex.pcic/issues
.. _contributor's guide: https://github.com/pacificclimate/climdex.pcic/blob/master/CONTRIBUTING.rst

Still need help?
----------------

* Contact climate@uvic.ca and let us know what we can do.

.. [1] Please know that the pool of people who can provide support for the package is extremely small and time is limited.  We don't necessarily have the capacity for long, open-ended user support. If you keep your questions short, specific and direct, there's a greater probability that someone will take on the ticket.
