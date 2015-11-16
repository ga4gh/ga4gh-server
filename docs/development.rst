.. _development:

-----------
Development
-----------

Thanks for your interest in helping us develop the GA4GH reference
implementation! There are lots of ways to contribute, and it's easy
to get up and running. This page should provide the basic information
required to get started; if you encounter any difficulties
`please let us know <https://github.com/ga4gh/server/issues>`_

.. warning::

    This guide is a work in progress, and is incomplete.

***********************
Development environment
***********************

We need a development Python 2.7 installation, Git, and some basic
libraries. On Debian or Ubuntu, we can install these using

.. code-block:: bash

    $ sudo apt-get install python-dev git zlib1g-dev libxslt1-dev libffi-dev libssl-dev

.. note::
    TODO: Document this basic step for other platforms? We definitely want
    to tell people how to do this with Brew or ports on a Mac.

If you don't have admin access to your machine, please contact your system
administrator, and ask them to install the development version of Python 2.7
and the development headers for `zlib <http://www.zlib.net/>`_.

Once these basic prerequisites are in place, we can then bootstrap our
local Python installation so that we have all of the packages we require
and we can keep them up to date. Because we use the functionality
of the recent versions of ``pip`` and other tools, it is important to
use our own version of it and not any older versions that may be
already on the system.

.. code-block:: bash

    $ wget https://bootstrap.pypa.io/get-pip.py
    $ python get-pip.py --user

This creates a `user specific <https://www.python.org/dev/peps/pep-0370/>`_
site-packages installation for Python, which is based in your ``~/.local``
directory. This means that you can now install any Python packages you like
without needing to either bother your sysadmin or worry about breaking your
system Python installation. To use this, you need to add the newly installed
version of ``pip`` to your ``PATH``. This can be done by adding something
like

.. code-block:: bash

    export PATH=$HOME/.local/bin:$PATH

to your ``~/.bashrc`` file. (This will be slightly different if you use
another shell like ``csh`` or ``zsh``.)

We then need to activate this configuration by logging out, and logging back in.
Then, test this by running:

.. code-block:: bash

    $ pip --version
    pip 6.0.8 from /home/username/.local/lib/python2.7/site-packages (python 2.7)

We are now ready to start developing!

***************
GitHub workflow
***************

First, go to https://github.com/ga4gh/server and click on the 'Fork'
button in the top right-hand corner. This will allow you to create
your own private fork of the server project where you can work.
See the `GitHub documentation <https://help.github.com/articles/fork-a-repo/>`_
for help on forking repositories.
Once you have created your own fork on GitHub, you'll need to clone a
local copy of this repo. This might look something like:

.. code-block:: bash

    $ git clone git@github.com:username/server.git

We can then install all of the packages that we need for developing the
GA4GH reference server:

.. code-block:: bash

    $ cd server
    $ pip install -r requirements.txt --user

This will take a little time as the libraries that we require are
fetched from PyPI and built.

It is also important to set up an
`upstream remote <https://help.github.com/articles/configuring-a-remote-for-a-fork/>`_
for your repo so that you can sync up with the changes that other people
are making:

.. code-block:: bash

    $ git remote add upstream https://github.com/ga4gh/server.git

All development is done against the ``develop`` branch.  The ``stable``
branch is meant to be kept stable since it is the branch releases are
based on -- don't touch it!  These are the two mainline branches.
Our branching model is loosely based on the one described
`here <http://nvie.com/posts/a-successful-git-branching-model/>`__.

All development should be done in a topic branch.  That is, a branch
that the developer creates him or herself.  These steps will create
a topic branch (replace ``TOPIC_BRANCH_NAME`` appropriately):

.. code-block:: bash

    $ git fetch --all
    $ git checkout develop
    $ git merge --ff-only upstream/develop
    $ git checkout -b TOPIC_BRANCH_NAME

Topic branch names should include the issue number (if there is a tracked
issue this change is addressing) and provide some hint as to what the
changes include.  For instance, a branch that addresses the (imaginary)
tracked issue with issue number #123 to add more widgets to the code
might be named ``123_more_widgets``.

At this point, you are ready to start adding, editing and deleting files.
Stage changes with ``git add``.  Afterwards, checkpoint your progress by
making commits:

.. code-block:: bash

    $ git commit -m 'Awesome changes'

(You can also pass the ``--amend`` flag to ``git commit`` if you want to
incorporate staged changes into the most recent commit.)

Once you have changes that you want to share with others, push your
topic branch to GitHub:

.. code-block:: bash

    $ git push origin TOPIC_BRANCH_NAME

Then create a pull request using the GitHub interface.  This pull request
should be against the ``develop`` branch (this should happen automatically).

At this point, other developers will weigh in on your changes and will
likely suggest modifications before the change can be merged into
``develop``.  When you get around to incorporating these suggestions,
it is likely that more commits will have been added to the ``develop``
branch.  Since you (almost) always want to be developing off of the
latest version of the code, you need to perform a rebase to incorporate
the most recent changes from ``develop`` into your branch.

.. warning::

    We recommend against using ``git pull``.  Use ``git fetch`` and ``git
    rebase`` to update your topic branch against mainline branches
    instead.  See the :ref:`Git Workflow Appendix <git-appendix>` for
    elaboration.

.. code-block:: bash

    $ git fetch --all
    $ git checkout develop
    $ git merge --ff-only upstream/develop
    $ git checkout TOPIC_BRANCH_NAME
    $ git rebase develop

At this point, several things could happen.  In the best case, the rebase
will complete without problems and you can continue developing.  In other
cases, the rebase will stop midway and report a merge conflict.  That is,
git has determined that it is impossible for it to determine how to
combine the changes from the new commits in the ``develop`` branch and
your changes in your topic branch and needs manual intervention to
proceed.  GitHub has some
`documentation <https://help.github.com/articles/resolving-merge-conflicts-after-a-git-rebase/>`_ on how to resolve rebase merge conflicts.

Once you have updated your branch to the point where you think that you
want to re-submit the code for other developers to consider, push the
branch again, this time using the force flag:

.. code-block:: bash

    $ git push --force origin TOPIC_BRANCH_NAME

If you had tried to push the topic branch without using the force flag,
it would have failed.  This is because non-force pushes only succeed when
you are only adding new commits to the tip of the existing remote branch.
When you want to do something other than that, such as insert commits
in the middle of the branch history (what ``git rebase`` does), or modify a
commit (what ``git commit --amend`` does) you need to blow away the remote
version of your branch and replace it with the local version.  This is
exactly what a force push does.

.. warning::

    Never use the force flag to push to the ``upstream`` repository.  Never use
    the force flag to push to the ``develop`` or ``stable`` branches.  Only use
    the force flag on your repository and on your topic branches.
    Otherwise you run the risk of screwing up the mainline branches, which
    will require manual intervention by a senior developer and manual
    changes by every downstream developer.  That is a recoverable
    situation, but also one that we would rather avoid.  (Note: a hint that
    this has happened is that one of the above listed merge commands that
    uses the ``--ff-only`` flag to merge a remote mainline branch into a
    local mainline branch fails.)

One task that you might be asked to do before your topic branch can be
merged is "squashing your commits."  We want the git history to be clean
and informative, and we do that by crafting one and only one commit
message per logical change.  In the normal course of development (unless
one is constantly committing with the ``--amend`` flag) many intermediate
commits can be created that should be squashed down to (usually) one before
it can be merged.  Do this with (assuming you are in your topic branch):

.. code-block:: bash

    $ git rebase -i `git merge-base develop HEAD`

This will launch an editor that will give you control over how you want
to structure your commits.  Usually you just want to "pick" the first
commit and "squash" all of the subsequent commits, and then ensure that
the final commit message is clean (best practice is to give a short
summary of the change on the first line, a blank line, and then a more
detailed description of the change following, with the issue number
-- if there is one -- in the detailed description).  More information
about the interactive rebase process can be found
`here <https://help.github.com/articles/about-git-rebase/>`__.
Once the commits are to your liking, you can push the branch to your
remote repository (which will require a force push if you reordered
or deleted commits that existed in the remote version of the branch).

(It usually is a good idea to squash commits before rebasing your topic
branch on top of a mainline branch.  See the elaboration in the :ref:`Git
Workflow Appendix <git-appendix>` on this topic.)

Once your pull request has been merged into ``develop``, you can close
the pull request and delete the remote branch in the GitHub interface.
Locally, run this command to delete the topic branch:

.. code-block:: bash

    $ git branch -D TOPIC_BRANCH_NAME

Only the tip of the iceberg of git and GitHub has been covered in this
section, and much more can be learned by browsing their documentation.
For instance, get help on the ``git commit`` command by running:

.. code-block:: bash

    $ git help commit

To master git, we recommend reading this free book (save chapter four,
which is about git server configuration): `Pro Git
<https://git-scm.com/book/en/v2>`_.


************
Contributing
************

See the files ``CONTRIBUTING.md`` and ``STYLE.md`` for an overview of
the processes for contributing code and the style guidelines that we
use.


*********************
Development utilities
*********************

All of the command line interface utilities have local scripts
that simplify development: for example, we can run the local version of the
``ga2sam`` program by using::

    $ python ga2sam_dev.py

To run the server locally in development mode, we can use the ``server_dev.py``
script, e.g.::

    $ python server_dev.py

will run a server using the default configuration. This default configuration
expects a data hierarchy to exist in the ``ga4gh-example-data`` directory.
This default configuration can be changed by providing a (fully qualified)
path to a configuration file (see the :ref:`configuration`
section for details).

There is also an OpenID Connect (oidc) provider you can run locally for
development and testing. It resides in ``/oidc-provider`` and has a run.sh
file that creates a virtualenv, installs the necessary packages, and
runs the server. Configuration files can be found in
``/oidc-provider/simple_op``::

    $ cd oidc-provider
    $ ./run.sh

The provider expects OIDC redirect URIs to be over HTTPS, so if the ga4gh
server is started with OIDC enabled, it defaults to HTTPS. You can run the
server against this using::

    $ python server_dev.py -c LocalOidConfig

************
Organisation
************

The code for the project is held in the ``ga4gh`` package, which corresponds to
the ``ga4gh`` directory in the project root. Within this package, the
functionality is split between the ``client``, ``server``, ``protocol`` and
``cli`` modules.  The ``cli`` module contains the definitions for the
``ga4gh_client`` and ``ga4gh_server`` programs.

An important file in the project is ``ga4gh/_protocol_definitions.py``.
This file defines the classes for the GA4GH protocol.
The file is generated using the ``scripts/process_schemas.py`` script,
which takes input data from the
`GA4GH schemas repo <https://github.com/ga4gh/schemas>`_.
To generate a new ``_protocol_definitions.py`` file, use

.. code-block:: bash

   $ python scripts/process_schemas.py -i path/to/schemas desiredVersion

Where ``desiredVersion`` is the version that will be written to the
``_protocol_definitions.py`` file.  This version must be in the form
``major.minor.revision`` where major, minor and revision can be any
alphanumeric string.

.. _git-appendix:

*********************
Git Workflow Appendix
*********************

++++++++++++++++++++++
Don't use ``git pull``
++++++++++++++++++++++

We recommend against using ``git pull``.  The ``git pull`` command by
default combines the ``git fetch`` and the ``git merge`` command.  If your
local branch has diverged from its remote tracking branch, running ``git
pull`` will create a merge commit locally to join the two branches.

In some workflows, this is not an issue.  For us, however, it creates a
problem in the future.  When you are ready to submit your topic branch in a
pull request, we ask you to squash your commits (usually down to one
commit).  Given the complex graph topography created by all of the merges, the
order in which git applies commits in the squash is very difficult to
reason about and will likely create merge conflicts that you find
unnecessary and nonsensical (and therefore, highly aggravating!).

We instead recommend using ``git fetch`` and ``git rebase`` to update your
local topic branch against a mainline branch.  This will create a linear
commit history in your topic branch, which will be easy to squash, since the
commits are applied in the squash in the order that you made them.

``git pull`` does have the ``--rebase`` option which will do a rebase
instead of a merge to incorporate the remote branch.  One can also set the
``branch.autosetuprebase always`` config option to have ``git pull`` do a
rebase by default (i.e. without passing the ``--rebase`` flag) instead of a
merge.  This will avoid the issue of squashing a non-linear commit history.

So, in truth, we are really recommending against squashing local branches
with many merge commits in them.  However, using the default settings for
``git pull`` is the easiest way to end up in this situation.  Therefore,
don't use ``git pull`` unless you know what you are doing.

+++++++++++++++++++
Squash, then rebase
+++++++++++++++++++

When updating a local topic branch with changes from a mainline branch, we
recommend squashing commits in your topic branch down to one commit before
rebasing on top of the mainline branch.  The reason for this is that, under the
hood, to apply the rebase ``git rebase`` essentially cherry-picks each
commit from your topic branch not in the mainline branch and applies it to the
mainline branch.  Each one of these applications can cause a merge
conflict.  It is much better to face the potential of only one merge
conflict than N merge conflicts (where N is the number of unique commits in the
local branch)!

The difficulty of proceeding the opposite way (rebasing, then squashing) is
only compounded because of the unintuitiveness of the N merge conflicts.
When presented with a merge conflict, your likely intuition is to put the
file in the state that you think it ought to be in, namely the condition it was
in after the Nth commit.  However, if that state was different than the
state that git thinks it should be in -- namely, the state of the file at
commit X where X<N -- then you have only created the potential for more
merge conflicts.  When the next intermediate commit, Y (where X<Y<N) is
applied, it too will create a merge conflict.  And so on.

So squash, then rebase, and avoid this whole dilemma.  The terms are a bit
confusing since both "squashing" and "rebasing" are accomplished via the
``git rebase`` command.  As mentioned above, squash the commits in your
topic branch with (assuming you have branched off of the ``develop``
mainline branch):

.. code-block:: bash

    $ git rebase -i `git merge-base develop HEAD`

(``git merge-base develop HEAD`` specifies the most recent commit that both
``develop`` and your topic branch share in common.  Normally this is
equivalent to the most recent commit of ``develop``, but that's not
guaranteed -- for instance, if you have updated your local ``develop``
branch with additional commits from the remote ``develop`` since you
created your topic branch which branched off of the local ``develop``.)

And rebase with (again, assuming ``develop`` as the mainline branch):

.. code-block:: bash

    $ git rebase develop

++++++++++++++++++++++++++++++
GitHub's broken merge/CI model
++++++++++++++++++++++++++++++

GitHub supports continuous integration (CI) via `Travis CI
<https://travis-ci.com/>`_.  On every pull request, Travis runs a suite of
tests to determine if the PR is safe to merge into the mainline branch that it
targets.  Unfortunately, the way that GitHub's merge model is structured
does not guarantee this property.  That is, it is possible for a PR to pass the
Travis tests but for the mainline branch to fail them after that PR is
merged.

How can this happen?  Let's illustrate by example: suppose PR A and PR B
both branch off of commit M, which is the most recent commit in the
mainline branch.  A and B both pass CI, so it appears that it is safe to
merge them into the mainline branch.  However, it is also true that the
changes in A and B have never been tested `together` until CI is run on the
mainline branch after both have been merged.  If PR A and B have
incompatible changes, even if both merge cleanly, CI will fail in the
mainline branch.

GitHub could solve this issue by not allowing a PR to be merged unless it
both passed CI and its branch contained (in addition to the commits it
wanted to merge in to mainline) every commit in the mainline branch.  That is,
no PR could be merged into mainline unless its commits were tested with
every commit already in mainline.  Right now GitHub does not mandate this
strict sequencing of commits, which is why it can never guarantee that the
mainline CI will pass, even if all the PR CIs passed.

Developers could also enforce this property manually, but we have
determined that not using GitHub's UI merging features and judiciously
re-submitting PRs for additional CI would be more effort than fixing a
broken test in a mainline branch once in a while.
