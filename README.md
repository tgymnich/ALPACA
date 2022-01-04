# *ALPACA*
```
\\
l '>
| |
| |
| ALPACA~
||    ||
''    ''
```
*ALPACA* is a MPI-parallelized C++ code framework to simulate compressible multiphase flow physics. It allows for advanced high-resolution sharp-interface modeling empowered with efficient multiresolution compression. The modular code structure offers a broad flexibility to select among many most-recent numerical methods covering WENO/T-ENO, Riemann solvers (complete/incomplete), strong-stability preserving Runge-Kutta time integration schemes, level-set methods and many more.  

⚠️ On macOS you have to increase the stack size by running `ulimit -s 65532`. This only mitigates the problem to a certain degree. If the new stack limit is reached ALPACA will still crash ⚠️

# Install

## macOS

```
brew install tgymnich/tap/hdf5-mpich
brew install cmake
brew install paraview 

git clone --recurse-submodules -j8 https://github.com/tgymnich/ALPACA.git
cd ALPACA
mkdir build && cd build
cmake ..
cmake --build .

ulimit -s 65532
mpirun -np 4 ./ALPACA ../inputfile.xml
```

# Terms of usage

*ALPACA* is free software (GNU GPLv3), see the [LICENSE file](LICENSE) and/or sourcefile header in the repository. If you use *ALPACA*, please cite the following paper:

```
@article{Hoppe2020,
title = "A modular massively parallel computing strategy for three-dimensional multiresolution simulations of compressible flows",
journal = "tba",
volume = "tba",
pages = "tba",
year = "2020",
doi = "https://doi.org/xxxx",
author = "Nils Hoppe and Stefan Adami and Nikolaus A. Adams",
keywords = "Multiresolution, Morton order, Compressible flow solver"
}
```

# First steps
To get started with *ALPACA*, please obtain a copy of the repository. 

## Using GIT:
If not available - install GIT on your machine.

Next, navigate to the place where you would like to have the ALPACA folder.

We include the third party libraries used within ALPACA via git submodules. Thus, ALPACA has to be cloned with all its submodules. This can be realized with the following git command: 
### Via `ssh`
``` 
git clone --recursive git@gitlab.lrz.de:nanoshock/ALPACA.git
```
### Via `https`
```
git clone --recursive https://gitlab.lrz.de/nanoshock/ALPACA.git
```

### Existing *ALPACA*
In case the submodules should be included in an already existing ALPACA git environment, all submodules can be updated using:
```
git submodule update --init --recursive
```

## Without GIT (not recommended):
Alternatively, you can directly download and unpack a compressed source-code package (zip / tar.gz / tar.bz2 / tar) from the main [project-page](./) in gitlab. Note, for this option you need to provide the third-party libraries separately.

## Compilation and first simulation
Now, follow the [installation instructions](../../wikis/1_Installation) to get your artiodactyla running. *ALPACA* is shipped with ready-to-use [examples](../../wikis/2_Examples).


# Documentation
The source code of *ALPACA* is fully documented, please also see the respective journal articles for further details.
Optionally, Doxygen can be used to generate a full documentation of the code.
The configuration file ["Doxyfile"](Doxyfile) is provided in the repository.
To generate the documentation (assuming "doxygen" is installed on your system), please use the following command:
```
> doxygen Doxyfile
```
Note, we used doxygen-1.8.11 to write this guide.

# Testing

*ALPACA* provides different levels of testing. Unit tests allow to test individual units of source code. Running the
unit tests for *ALPACA* is fast and takes only about a minute. We also provide a Testsuite for more extensive tests.
In the following, details on how to run the respective tests are given.

## Unit Tests
To compile and run the unit tests the following steps are necessary:
* Go to the *ALPACA* repository folder
```
cd /path/to/ALPACA
```
* Create a build folder and run cmake
```
mkdir build
cd build
cmake ..
```
* Create the unit test executable *Paco*
```
make Paco -j 4
```
* Run the single- as well as two-core unit tests
```
mpiexec -n 1 ./Paco [1rank]
mpiexec -n 2 ./Paco [2rank]
```

## Testsuite
The Testsuite serves to validate successful completion of several single- and two-phase cases, MPI parallelization, symmetry
preserving techniques as well as correct physical behaviour. In order to run the Testsuite a python interpreter
(Python >= 3.6) is required. Then, the following steps are necessary:

* Go to the *ALPACA* repository folder
```
cd /path/to/ALPACA
```
* Go to the folder where the `run_testsuite.py` file is located
```
cd Utilities/TestSuite
```
* Run the Testsuite with a selected configuration
```
python run_testsuite.py ../Inputfiles/TestSuiteConfiguration/AerConfig.xml
```
The Testsuite takes several hours to complete and requires non-negligible computational resources (several cores).
For testing during development we recommend to adjust the extent of the considered tests via the configuration file 
(see `Utilities/Inputfiles/TestSuiteConfiguration/AerConfig.xml`).

Running the Testsuite on the [Linux cluster of the LRZ](https://doku.lrz.de/display/PUBLIC/Linux+Cluster) is done in two steps. The first step is to compile all
relevant executables:
* Go to the *ALPACA* repository folder
```
cd /path/to/ALPACA
```
* Go to the Testsuite folder
```
cd Utilities/TestSuite
```
* Run `compile_merge_suite.sh` (`nohup` can be used to detach the process from the terminal)
```
bash compile_merge_suite.sh
```
After this step is completed (takes several hours) submit a job to run the Testsuite:
* Go to the *ALPACA* repository folder
```
cd /path/to/ALPACA
```
* Go to the Testsuite folder
```
cd Utilities/TestSuite
```
* Adapt the `testsuite.job` Slurm job file according to your needs and submit the job
```
sbatch testsuite.job
```
After the job has finished (again takes several hours) check the output file of the Slurm job for successful completion of the Testsuite.

# Wiki

The [wiki pages](../../wikis) contain introductory information on [installation instructions](../../wikis/1_Installation) and first [examples](../../wikis/2_Examples) to try out *ALPACA*. There is also some more background information on the [numerics](../../wikis/3_Numerics) together with a [developer's corner](../../wikis/4_DevelopersCorner) for experts.


# Collaborations

We are highly interested in fruitful collaborations and hope to provide a
useful tool to other research groups and interested scientists. If you are
working with *ALPACA*, we highly appreciate your comments and experiences to
improve the code. If you work on new features, please feel free to 
[contact us](mailto:nanoshock@aer.mw.tum.de?subject=Interest%20on%20collaboration)
to avoid redundant developments.

# Q&A

If you encounter any problems with *ALPACA* regarding e.g. compilation or
performing your simulations, please don't hesitate to contact the developers
via

1.  Gitlab's 'Issue tracking system' (recommended) or
2.  get in touch with us by [mail](mailto:nanoshock@aer.mw.tum.de?subject=Question%20regarding%20ALPACA)

# Peer-reviewed publications using Alpaca

> [Winter, J. M., Kaiser, J. W. J., Adami, S., & Adams, N. A. (2019). Numerical investigation of 3D drop-breakup mechanisms using a sharp interface level-set method. Proceedings of the 11th Symposium on Turbulence and Shear Flow Phenomena, Southampton, United Kingdom.](https://mediatum.ub.tum.de/1522845)

> [Kaiser, J. W. J., Adami, S., & Adams, N. A. (2019). Three-dimensional direct numerical simulation of shock-induced bubble collapse near gelatin. Proceedings of the 11th Symposium on Turbulence and Shear Flow Phenomena, Southampton, United Kingdom.](https://mediatum.ub.tum.de/1522844)

> [Kaiser, J. W. J., Hoppe, N., Adami, S., & Adams, N. A. (2019). An adaptive local time-stepping scheme for multiresolution simulations of hyperbolic conservation laws. Journal of Computational Physics: X, 4, 100038.](https://www.sciencedirect.com/science/article/pii/S259005521930054X)

<!---
> [Paula, T., Adami, S. & Adams, N. A. (2019). Analysis of the early stages of liquid-water-drop explosion by numerical simulation. Phys. Rev. Fluids, 4(4), 044003.](https://link.aps.org/doi/10.1103/PhysRevFluids.4.044003)
-->

> [Fleischmann, N., Adami, S. & Adams, N. A. (2019). Numerical symmetry-preserving techniques for low-dissipation shock-capturing schemes. Computers & Fluids, 189, 94-107.](https://doi.org/10.1016/j.compfluid.2019.04.004)

> [Fleischmann, N., Adami, S. & Adams, N. A. (2020). A low dissipation method to cure the grid-aligned shock instability. Journal of Computational Physics, 401, 109004.](https://mediatum.ub.tum.de/1536321)

> [Kaiser, J. W. J., Adami, S., Akhatov, I. S. & Adams, N. A. (2020). A semi-implicit conservative sharp-interface method for liquid-solid phase transition. International Journal of Heat and Mass Transfer, 155, 119800.](https://doi.org/10.1016/j.ijheatmasstransfer.2020.119800)


## Bibtex keys for publications
```
@inproceedings{Winter2017,
author = " Winter, J.M., Kaiser, J.W.J., Adami, S. and Adams, N.A.", 
title = "Numerical investigation of 3D drop-breakup mechanisms using a sharp interface level-set method",
booktitle = "11th International Symposium on Turbulence and Shear Flow Phenomena, TSFP 2019",
address = "Southampton; United Kingdom",
year = "2019",
url = "https://mediatum.ub.tum.de/1522845"
}

@inproceedings{Kaiser2019 ,
author = " Kaiser, J.W.J., Adami, S. and Adams, N.A.", 
title = "Three-dimensional direct numerical simulation of shock-induced bubble collapse near gelatin",
booktitle = "11th International Symposium on Turbulence and Shear Flow Phenomena, TSFP 2019",
address = "Southampton; United Kingdom",
year = "2019",
url ="https://mediatum.ub.tum.de/1522844"
}

@article{Kaiser2019b,
author = " Kaiser, J.W.J., Hoppe, N., Adami, S. and Adams, N.A.", 
title = "An adaptive local time-stepping scheme for multiresolution simulations of hyperbolic conservation laws",
journal = "Journal of Computational Physics: X",
volume = "4",
pages = "100038",
year = "2019",
issn = "2590-0552",
doi = "https://doi.org/10.1016/j.jcpx.2019.100038",
url = "http://www.sciencedirect.com/science/article/pii/S259005521930054X"
}

@article{Fleischman2019,
author = "Fleischmann, N., Adami, S. and Adams, N.A.",
title = "Numerical symmetry-preserving techniques for low-dissipation shock-capturing schemes",
journal = "Computers & Fluids",
volume = "189",
pages = "94 - 107",
year = "2019",
issn = "0045-7930",
doi = "https://doi.org/10.1016/j.compfluid.2019.04.004",
url = "http://www.sciencedirect.com/science/article/pii/S0045793018308399",
}

@article{Fleischmann2020,
author = "Fleischmann, N., Adami, S. and Adams, N.A.",
title = "A low dissipation method to cure the grid-aligned shock instability",
journal = "Journal of Computational Physics",
volume = "401",
pages = "109004",
year = "2020",
issn = "0021-9991",
doi = "https://doi.org/10.1016/j.jcp.2019.109004"
}

@article{Kaiser2020,
author = "Kaiser, J.W.J., Adami, S., Akhatov, I.S. and Adams, N.A.",
title = "A semi-implicit conservative sharp-interface method for liquid-solid phase transition",
journal = "International Journal of Heat and Mass Transfer",
volume = "155",
pages = "119800",
year = "2020",
issn = "0017-9310",
doi = "https://doi.org/10.1016/j.ijheatmasstransfer.2020.119800",
url = "http://www.sciencedirect.com/science/article/pii/S0017931019361873"
}
```

# Acknowledgments

*ALPACA* was fed and grown thanks to several supporter:
* This project has received funding from the European Research Council (ERC) under the European Union’s Horizon 2020 research and innovation program: ERC Advanced Grant No. 667483, Prof. Dr. Nikolaus A. Adams, "NANOSHOCK - Manufacturing Shock Interactions for Innovative Nanoscale Processes"
* This project has received computing time on the GCS Supercomputer SuperMUC at Leibniz Supercomputing Centre (www.lrz.de) from the Gauss Centre for Supercomputing e.V. (www.gauss-centre.eu).
* This project has received funding from the Bavarian State Ministry of Science and the Arts through the Competence Network for Scientific High Performance Computing in Bavaria (KONWIHR).
* This project has received funding from German Research Foundation (DFG).
