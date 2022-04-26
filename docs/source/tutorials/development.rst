===========
Development
===========

How to contribute
=================

Development of new features should be done in separate branches (even on forks)
of the repository in order to ensure that:

- the master branch never contains unfinished/unvalidated code
- unnecessary build notifications to the slack channel are avoided

PULL REQUESTS. (Code review)

Automated testing
=================

TRAVIS

Code quality
============

While measures of code quality (and aesthetic) are highly subjective, a common
standard helps for the collaborative development of open-source code. As most
Python projects molSimplify follows the recommendations outlined in PEPXXX.

Linter
------

The easiest way to ensure your contributions conform to PEPXXX is to install
a linter such as flake8 for your IDE of choice.

EXPLAIN what the github action checks for.

Pre-commit hooks
----------------

In order to automate the process of code lintering molSimplify also provides a
configuration file for a pre-commit hook.