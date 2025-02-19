
auto-AUTO (or AUTO²)
====================

General Information
-------------------

AUTO² or auto-AUTO is an [AUTO](https://github.com/auto-07p/auto-07p) automatic search algorithm codebase
to enhance the original AUTO-07p Python interface with a top layer which allows users to:

* automate the continuation of as many branches as possible, branching whenever possible to construct full
  bifurcation trees, and finishing computations based on a predefined logic
  (meeting other branches, looping branches, etc...)
* plot results with |Matplotlib|
* perform these computations in [Jupyter](https://jupyter.org) notebooks

About
-----

(c) 2025 Jonathan Demaeyer and Oísin Hamilton. 

See [LICENSE.txt](https://raw.githubusercontent.com/Climdyn/auto-AUTO/master/LICENSE.txt) for license information.

Please consult the auto-AUTO [code repository](http://www.github.com/Climdyn/auto-AUTO) for updates.

Installation
------------

#### Installing AUTO

> To use auto-AUTO, you need the [bleeding edge version of AUTO](https://github.com/auto-07p/auto-07p) available 
> on Github for this codebase to work properly !

Here how to install AUTO from GitHub:

* First clone the AUTO repository somewhere:


    git clone https://github.com/auto-07p/auto-07p.git

* Then in a terminal, in the created folder, run:


    ./configure
    make

* Your AUTO installation should now be finished, but you still need to add the 
following line to your `.bashrc` file:


    source [path-to-auto-07p]/cmds/auto.env.sh

In addition, we recommend that you edit this file so that the `AUTO_DIR` environment 
variable specified there points to the correct folder where you installed AUTO.

 > Be sure to have all the AUTO requirements pre-installed. See AUTO documentation for 
> more details. In case of issues, we recommend reading the documentation completely.

After that last step, you have AUTO properly configured and are ready to install auto-AUTO.

> If AUTO version is changing over time, you need to update the version from GitHub and do
> the installation again.

#### Installing auto-AUTO with pip

TODO

#### Installing auto-AUTO with Anaconda

The second-easiest way to install and run qgs is to use an appropriate 
environment created through [Anaconda](https://www.anaconda.com/).

First install Anaconda and clone the repository:

    git clone https://github.com/Climdyn/auto-AUTO.git

Then install and activate the Python3 Anaconda environment:

    conda env create -f environment.yml
    conda activate auto2

and the code is installed. You can test the Jupyter notebooks present in the 
[notebooks folder](./notebooks).
For instance, running

    conda activate auto2
    cd notebooks
    jupyter-notebook

will lead you to your favorite browser where you can load and run the examples.

Documentation
-------------

To build the documentation, please run (with the conda environment activated):

    cd documentation
    make html

You may need to install [make](https://www.gnu.org/software/make/) if it is not already present on your system.
Once built, the documentation is available [here](./documentation/build/html/index.html).

Forthcoming developments
------------------------

* Regime diagrams object

Contributing to qgs
-------------------

If you want to contribute actively, please contact the main authors.

In addition, if you have made changes that you think will be useful to others, please feel free to suggest these as a pull request on the [auto-AUTO Github repository](https://github.com/Climdyn/auto-AUTO).
