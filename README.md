# homeR

`homeR` is a collection of functions that we have
found useful at Neurobat in the analysis of experiments on buildings.

## Installation

Get the released version from CRAN:

```r
install.packages("homeR")
```

Or the development version from GitHub:

```r
# install.packages("devtools")
devtools::install_github("neurobat/homeR")
```

## Building the package

We recommend using RStudio to build this package. Make sure you have the latest
version installed (I'm currently using 0.99.896), and make also sure you have
`roxygen2` and `testthat` installed:

```r
install.packages("roxygen2")
install.packages("testthat")
```

From RStudio, create a new project (File -> New Project -> Existing Directory) and give
the project directory. Configure RStudio to build the documentation from Roxygen:
Build -> Configure Build Tools -> Generate documentation with Roxygen. Under
"Configure" ensure that "Rd files" and "NAMESPACE file" are checked.

Verify that you can build the package by running Document (Shift-Ctrl/Cmd-D) and
Test Package (Shift-Ctrl/Cmd-T).

## Committing generated files

The `NAMESPACE` and `*.Rd` files, although generated from the sources, should be
committed and kept under Git so that anyone cheking out the project from GitHub
can build and install it from the command line.