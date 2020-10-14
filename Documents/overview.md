# SABER overview
## Directory structure
The SABER repository is organized as follows:
* **CI**: continuous integration setup files and scripts
* **cmake**: SABER-specific compiler flags
* **Documents**: documentation
  * **autodoc**: code auto-documentation
* **src**: source code
  * **saber**: SABER source code
    * **bump**: BUMP core and interfaces
    * **external**: external tools
    * **gaugrid**: Gaussian grid tools
    * **interpolation**: interpolation interface
    * **oops**: OOPS interface
    * **util**: shared tools
* **test**: testing tools and source code
  * **mains**: executables sources
    * **model**: model-specific interfaces for the BUMP executable
  * **testinput**: input YAML files
  * **testlist**: lists of tests, data and reference files
* **tools**: useful bash and python scripts


## Code size and characteristics
Code report obtained with [CLOC](https://github.com/AlDanial/cloc).

### Documents

| language | files | blank | comment | code | comment/code ratio |
|:--------:|:--------:|:--------:|:--------:|:--------:|:--------:|
| Markdown | 51 | 52 | 0 | 772 | 0 % |
| C/C++ Header | 1 | 2 | 21 | 3 | 700 % |


### src/saber/bump

| language | files | blank | comment | code | comment/code ratio |
|:--------:|:--------:|:--------:|:--------:|:--------:|:--------:|
| Fortran 90 | 35 | 5153 | 8066 | 20475 | 39 % |
| C/C++ Header | 1 | 6 | 7 | 39 | 17 % |
| C | 1 | 2 | 6 | 21 | 28 % |


### src/saber/external

| language | files | blank | comment | code | comment/code ratio |
|:--------:|:--------:|:--------:|:--------:|:--------:|:--------:|
| Fortran 90 | 4 | 366 | 1672 | 1766 | 94 % |


### src/saber/gaugrid

| language | files | blank | comment | code | comment/code ratio |
|:--------:|:--------:|:--------:|:--------:|:--------:|:--------:|
| Fortran 90 | 1 | 59 | 87 | 102 | 85 % |


### src/saber/interpolation

| language | files | blank | comment | code | comment/code ratio |
|:--------:|:--------:|:--------:|:--------:|:--------:|:--------:|
| Fortran 90 | 2 | 165 | 231 | 407 | 56 % |
| C/C++ Header | 2 | 23 | 19 | 56 | 33 % |
| C++ | 1 | 12 | 17 | 41 | 41 % |


### src/saber/oops

| language | files | blank | comment | code | comment/code ratio |
|:--------:|:--------:|:--------:|:--------:|:--------:|:--------:|
| C/C++ Header | 13 | 343 | 235 | 1220 | 19 % |


### src/saber/util

| language | files | blank | comment | code | comment/code ratio |
|:--------:|:--------:|:--------:|:--------:|:--------:|:--------:|
| Fortran 90 | 9 | 1170 | 1641 | 2692 | 60 % |
| Fortran 77 | 2 | 28 | 37 | 102 | 36 % |


### test/mains

| language | files | blank | comment | code | comment/code ratio |
|:--------:|:--------:|:--------:|:--------:|:--------:|:--------:|
| Fortran 90 | 2 | 168 | 279 | 725 | 38 % |
| C++ | 5 | 12 | 41 | 67 | 61 % |


## Code auto-documentation
### Documents: documentation

| Name | Purpose |
| :--: | :---------- |


### src/saber/bump: BUMP core and interfaces

| Name | Purpose |
| :--: | :---------- |
| [tools_fit](autodoc/tools_fit.md) | fit-related tools |
| [tools_func](autodoc/tools_func.md) | usual functions |
| [tools_samp](autodoc/tools_samp.md) | sampling functions |
| [type_avg_blk](autodoc/type_avg_blk.md) | averaged statistics block derived type |
| [type_avg](autodoc/type_avg.md) | average routines |
| [type_bpar](autodoc/type_bpar.md) | block parameters derived type |
| [type_bump](autodoc/type_bump.md) | BUMP derived type |
| [type_bump_interface](autodoc/type_bump_interface.md) | BUMP derived type interface |
| [type_cmat_blk](autodoc/type_cmat_blk.md) | correlation matrix derived type |
| [type_cmat](autodoc/type_cmat.md) | C matrix derived type |
| [type_com](autodoc/type_com.md) | communications derived type |
| [type_cv_blk](autodoc/type_cv_blk.md) | control vector derived type |
| [type_cv](autodoc/type_cv.md) | control vector derived type |
| [type_diag_blk](autodoc/type_diag_blk.md) | diagnostic block derived type |
| [type_diag](autodoc/type_diag.md) | diagnostic derived type |
| [type_ens](autodoc/type_ens.md) | ensemble derived type |
| [type_geom](autodoc/type_geom.md) | geometry derived type |
| [type_hdiag](autodoc/type_hdiag.md) | hybrid diagnostics derived type |
| [type_io](autodoc/type_io.md) | I/O derived type |
| [type_lct_blk](autodoc/type_lct_blk.md) | LCT data derived type |
| [type_lct](autodoc/type_lct.md) | LCT data derived type |
| [type_linop](autodoc/type_linop.md) | linear operator derived type |
| [type_mesh](autodoc/type_mesh.md) | mesh derived type |
| [type_minim](autodoc/type_minim.md) | minimization data derived type |
| [type_mom_blk](autodoc/type_mom_blk.md) | moments block derived type |
| [type_mom](autodoc/type_mom.md) | moments derived type |
| [type_nam](autodoc/type_nam.md) | namelist derived type |
| [type_nicas_blk](autodoc/type_nicas_blk.md) | NICAS data block derived type |
| [type_nicas](autodoc/type_nicas.md) | NICAS data derived type |
| [type_obsop](autodoc/type_obsop.md) | observation operator data derived type |
| [type_samp](autodoc/type_samp.md) | sampling derived type |
| [type_tree](autodoc/type_tree.md) | tree derived type |
| [type_var](autodoc/type_var.md) | variance derived type |
| [type_vbal_blk](autodoc/type_vbal_blk.md) | vertical balance block derived type |
| [type_vbal](autodoc/type_vbal.md) | vertical balance derived type |


### src/saber/external: external tools

| Name | Purpose |
| :--: | :---------- |
| [tools_asa007](autodoc/tools_asa007.md) | inverse of symmetric positive definite matrix routines |
| [tools_qsort](autodoc/tools_qsort.md) | qsort routines |
| [tools_sp](autodoc/tools_sp.md) | spectral transforms |
| [tools_stripack](autodoc/tools_stripack.md) | STRIPACK routines |


### src/saber/gaugrid: Gaussian grid tools

| Name | Purpose |
| :--: | :---------- |
| [type_gaugrid](autodoc/type_gaugrid.md) | Gaussian grid type |


### src/saber/interpolation: interpolation interface

| Name | Purpose |
| :--: | :---------- |


### src/saber/oops: OOPS interface

| Name | Purpose |
| :--: | :---------- |


### src/saber/util: shared tools

| Name | Purpose |
| :--: | :---------- |
| [tools_atlas](autodoc/tools_atlas.md) | random numbers generator derived type |
| [tools_const](autodoc/tools_const.md) | define usual constants and missing values |
| [tools_kinds](autodoc/tools_kinds.md) | kinds definition |
| [tools_repro](autodoc/tools_repro.md) | reproducibility functions |
| [tools_atlas](autodoc/tools_atlas.md) | random numbers generator derived type |
| [type_mpl](autodoc/type_mpl.md) | MPI parameters derived type |
| [type_msv](autodoc/type_msv.md) | deal with missing values |
| [type_rng](autodoc/type_rng.md) | random numbers generator derived type |
| [type_timer](autodoc/type_timer.md) | timer data derived type |


### test/mains: executables sources

| Name | Purpose |
| :--: | :---------- |
| [type_model](autodoc/type_model.md) | model routines |


