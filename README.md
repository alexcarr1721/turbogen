# turbogen: Generate turbulent field on Cartesian grid

> Generates an isotropic turbulent field in a 3 dimensional rectangular grid.

---

## Table of Contents

- [Installation](#installation)
- [Overview](#overview)
- [Features](#features)
- [Support](#support)
- [License](#license)


---

## Installation

### Linux

- Clone this repo to a directory on your local machine that you have read/write priviliges using the following command,

```shell
$ git clone https://github.com/alexcarr1721/turbogen
```

- Make sure that you have an up to date version of cmake installed, along with MPI compilers, a parallel HDF5 build, and the intel MKL library. If you do not have all of these installed on your system then please do so before moving on to the next step. Make sure to add the binary directories for each of the above in your path, as well as the library directories in your $LD_LIBRARY_PATH variable. To build the turbogen program, do the following from the top-level directory <installpath>/turbogen

```shell
$ cd build
$ export FC=<your_favorite_Fortran_compiler>
$ cmake .. <options>
$ make
```

---

## Overview

- See the documentation under /doc for more info.

## Features
## Usage (Optional)
## Documentation (Optional)
## Tests (Optional)

- Going into more detail on code and technologies used
- I utilized this nifty <a href="https://github.com/adam-p/markdown-here/wiki/Markdown-Cheatsheet" target="_blank">Markdown Cheatsheet</a> for this sample `README`.
---

## Support

Please email me with any questions. (alexcarr.1721@gmail.com)

---
