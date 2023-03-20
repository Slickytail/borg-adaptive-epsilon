## BORG MOEA: Adaptive Epsilon

COCO experiment for the Borg MOEA with automatic epsilon adaptation.
After compilation of coco, the files `borg.c`, `borg.h`, `mt19937ar.c`, `mt19937ar.h`, `coco.c`, and `coco.h` should be placed in the root directory.

Note: For the moment, it is necessary to modify `borg.h` to include the definitions of `struct BORG_Algorithm_t` and `struct BORG_Population_t`. This will be changed in a future release.
