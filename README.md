# Statlab
R library for Statistical Modelling

## Usage
Import this file in RStudio with and then call setup like this:
```
source('path_to_statlab.R')
sl_setup()
```
## Extending Statlab
Feel free to create new functions, open a pull request once you're sure they work and I'll add them. Just a few of caveats:
- Do not set seeds inside functions, this should be done explicitly by the caller.
- Names should start with sl_ to avoid ambiguity with already existing functions.
- Remember to document your functions and add a cat at the end of the library with a brief description to be printed at setup.
