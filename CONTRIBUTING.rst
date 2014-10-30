Contributing to the climdex.pcic R package
==========================================

Getting Started
---------------

- Create a `Github account`_.
- Fork the repository on Github at https://github.com/pacificclimate/climdex.pcic.
- Work on the code (see the next section)
- Send us a `pull request`_.

.. _Github account: https://github.com/signup/free
.. _pull request: https://help.github.com/articles/using-pull-requests/

How to set up a development environment
---------------------------------------

You don't need much to get started for development. You'll need to have installed:

- R (ensure that all of the "Depends", "Imports", and "Suggests" pacakges are also installed)
- Any C++ build environment supported by the `CRAN package checking`_
- git
- your text editor of choice

That's it!

Once you have the required software installed, create a local clone of the repository.
::
    $ git clone https://github.com/[your_user]/climdex.pcic.git

Make sure that everything builds out of the box
::
    $ R CMD build climdex.pcic/

.. _CRAN package checking: http://cran.r-project.org/web/checks/check_flavors.html

How to run the tests
--------------------

Running the tests can be done with one command:
::
    hiebert@medusa:/home/data/projects/R_packages$ R CMD check climdex.pcic

You'll see a bunch of package building spew that has nothing to do with the tests. But towards the end, you see something like this:
::
    * checking tests ...
    Running ‘bootstrap.R’
    Running ‘test_indices.R’
    Running ‘test_primitives.R’
    OK

Bug reports
-----------

If there are problems with our package or bugs in the code, please let us know! We welcome bug reports. To submit one:

- `Create a new issue`_ on our GitHub page.
- Tag/label the issue as a bug
- Leave it unassigned

Then please follow these guidelines for writing your report:

- Please describe in as much detail as possible
- Include a complete description of:

  - Exactly what you did (i.e. "steps to reproduce")
  - What you expected to happen?
  - What did happen?

- Include *all* output from the terminal.
- Run R's ``sessionInfo()`` function and include the full output.

I cannot stress enough how important it is to contrast what you expected to happen, with what actually happened. When executing the code does not produce the *advertised* result, there is a bug in the package. When the code does not produce the result that you *wished* it had, this is *not* a bug. We receive far too many reports in the latter category.

.. _Create a new issue: https://github.com/pacificclimate/climdex.pcic/issues/new

How to build the docs
---------------------

The package documentation is inline in the code. All of the manual pages are built by using ``roxygen2``. Make sure that you have ``roxygen2`` installed and loaded:
::
    hiebert@medusa:/home/data/projects/R_packages/climdex.pcic$ R

    R version 3.0.2 (2013-09-25) -- "Frisbee Sailing"
    Copyright (C) 2013 The R Foundation for Statistical Computing
    Platform: x86_64-pc-linux-gnu (64-bit)

    R is free software and comes with ABSOLUTELY NO WARRANTY.
    You are welcome to redistribute it under certain conditions.
    Type 'license()' or 'licence()' for distribution details.

    Natural language support but running in an English locale

    R is a collaborative project with many contributors.
    Type 'contributors()' for more information and
    'citation()' on how to cite R or R packages in publications.

    Type 'demo()' for some demos, 'help()' for on-line help, or
    'help.start()' for an HTML browser interface to help.
    Type 'q()' to quit R.

    > library(roxygen2)

Then call ``roxygenize()`` to build the docs.
::
    > roxygen2::roxygenize(clean=T)
    NOTE: arguments in definition for validity method for class 'climdexInput' changed from (x) to (object)
    Writing NAMESPACE
    Writing climdex.pcic.Rd
    Writing ec.1018935.Rd
    Writing get.last.monthday.of.year.Rd
    Writing climdexInput.Rd
    Writing get.outofbase.quantiles.Rd
    Writing climdexInput.raw.Rd
    Writing climdexInput.csv.Rd
    Writing climdex.fd.Rd
    Writing climdex.su.Rd
    Writing climdex.id.Rd
    Writing climdex.tr.Rd
    Writing climdex.gsl.Rd
    Writing climdex.txx.Rd

Submitting pull requests
------------------------

We would love help from the greater climate community in developing the package and we welcome contributions to climdex.pcic package.

- Please write tests for any functionality that you may add.
- Please modify tests for any functionality that you change.
- In short, please make sure that all of the tests pass.

After you are *positive* that everything is completely tested with passing test suite, we would love to see your pull request. If you are not familiar with the process, please follow the GitHub's help page for submitting `pull request`_.

Don't code? No problem!
-----------------------

Even if you don't program for a living there are plenty of ways to help. Not only is the code open and collaborative, but so is the documentation and issue tracking. Anyone can help with these. If you can't program, consider helping with the following:

- If the documentation doesn't answer your questions, it probably doesn't answer many people's questions. Help us all out and write something that does.
- Take a look through the outstanding `"help wanted" issues`_, and see if you know any of the answers.
- If there are `open bug reports`_, see if you can reproduce the problem and verify that it exists. Having bug reports validated and/or clarified by multiple parties is extremely valuable.
- Tell us your story. If ``climdex.pcic`` has helped your project to better understand climate extremes, we would love to hear about it. Write a blog post and/or send an e-mail to the `package maintainer`_.

.. _"help wanted" issues: https://github.com/pacificclimate/climdex.pcic/labels/help%20wanted
.. _open bug reports: https://github.com/pacificclimate/climdex.pcic/labels/bug
.. _package maintainer: mailto:hiebert@uvic.ca
