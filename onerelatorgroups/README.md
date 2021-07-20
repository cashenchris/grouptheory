# onerelatorgroups
Some scripts for computations in finitely generated one-relator groups.
First get:
https://github.com/cashenchris/grouptheory.git
https://github.com/cashenchris/freegroups.git
https://github.com/cashenchris/automaticgroups.git
https://github.com/cashenchris/smallcancellation.git

The main feature is the function certify_hyperbolicity to check if a one-relator group is hyperbolic or not. This checks several criteria. If the initial checks are all inconclusive then the function attempts to check using other programs, first using the "walrus" package of the GAP program, and second using the program kbmag. To perform these additional checks those programs must be installed and on the path.

GAP is available from:
https://www.gap-system.org/

kbmag is available from:
https://homepages.warwick.ac.uk/~mareg/download/kbmag2/

