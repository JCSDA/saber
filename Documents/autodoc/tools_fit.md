# Module tools_fit

| Type | Name | Purpose | Arguments          |
| :--: | :--: | :------ | :----------------- |
| subroutine | [fast_fit](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/tools_fit.F90#L37) | fast fit length-scale estimation based on the value at mid-height | <b>mpl</b> :  MPI data - type(mpl_type) - inout<br><b>n</b> :  Vector size - integer - in<br><b>iz</b> :  Zero separation index - integer - in<br><b>dist</b> :  Distance - real(kind_real) - in<br><b>raw</b> :  Raw data - real(kind_real) - in<br><b>fit_r</b> :  Fast fit result - real(kind_real) - out |
| subroutine | [ver_smooth](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/tools_fit.F90#L183) | homogeneous smoothing of a vertical profile | <b>mpl</b> :  MPI data - type(mpl_type) - inout<br><b>n</b> :  Vector size - integer - in<br><b>x</b> :  Coordinate - real(kind_real) - in<br><b>rv</b> :  Filtering support radius - real(kind_real) - in<br><b>profile</b> :  Vertical profile - real(kind_real) - inout |
| subroutine | [ver_fill](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/tools_fit.F90#L236) | missing values filling of a vertical profile | <b>mpl</b> :  MPI data - type(mpl_type) - inout<br><b>n</b> :  Vector size - integer - in<br><b>x</b> :  Coordinate - real(kind_real) - in<br><b>profile</b> :  Vertical profile - real(kind_real) - inout |
