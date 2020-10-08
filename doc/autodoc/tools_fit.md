# Module tools_fit

| Type | Name | Purpose | Arguments          |
| :--: | :--: | :------ | :----------------- |
| subroutine | [fast_fit](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/tools_fit.F90#L37) | fast fit length-scale estimation based on the value at mid-height | **mpl** :  MPI data - type(mpl_type) - inout<br>**n** :  Vector size - integer - in<br>**iz** :  Zero separation index - integer - in<br>**dist** :  Distance - real(kind_real) - in<br>**raw** :  Raw data - real(kind_real) - in<br>**fit_r** :  Fast fit result - real(kind_real) - out |
| subroutine | [ver_smooth](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/tools_fit.F90#L183) | homogeneous smoothing of a vertical profile | **mpl** :  MPI data - type(mpl_type) - inout<br>**n** :  Vector size - integer - in<br>**x** :  Coordinate - real(kind_real) - in<br>**rv** :  Filtering support radius - real(kind_real) - in<br>**profile** :  Vertical profile - real(kind_real) - inout |
| subroutine | [ver_fill](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/tools_fit.F90#L236) | missing values filling of a vertical profile | **mpl** :  MPI data - type(mpl_type) - inout<br>**n** :  Vector size - integer - in<br>**x** :  Coordinate - real(kind_real) - in<br>**profile** :  Vertical profile - real(kind_real) - inout |
