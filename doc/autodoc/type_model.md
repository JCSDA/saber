# Module type_model

| Type | Name | Purpose | Arguments          |
| :--: | :--: | :------ | :----------------- |
| subroutine | [model_alloc](https://github.com/JCSDA/saber/tree/develop/test/mains/type_model.F90#L142) | allocation | **model** :  Model - class(model_type) - inout |
| subroutine | [model_dealloc](https://github.com/JCSDA/saber/tree/develop/test/mains/type_model.F90#L165) | release memory | **model** :  Model - class(model_type) - inout |
| subroutine | [model_setup](https://github.com/JCSDA/saber/tree/develop/test/mains/type_model.F90#L212) | setup model | **model** :  Model - class(model_type) - inout<br>**mpl** :  MPI data - type(mpl_type) - inout<br>**nam** :  Namelist variables - type(nam_type) - inout |
| subroutine | [model_read](https://github.com/JCSDA/saber/tree/develop/test/mains/type_model.F90#L720) | read member field | **model** :  Model - class(model_type) - inout<br>**mpl** :  MPI data - type(mpl_type) - inout<br>**nam** :  Namelist - type(nam_type) - in<br>**filename** :  File name - character(len=*) - in<br>**fieldset** :  Fieldset - type(fieldset_type) - inout |
| subroutine | [model_read_member](https://github.com/JCSDA/saber/tree/develop/test/mains/type_model.F90#L771) | read member field | **model** :  Model - class(model_type) - inout<br>**mpl** :  MPI data - type(mpl_type) - inout<br>**nam** :  Namelist - type(nam_type) - in<br>**filename** :  File name - character(len=*) - in<br>**ie** :  Ensemble member index - integer - in<br>**fieldset** :  Fieldset - type(fieldset_type) - out |
| subroutine | [model_load_ens](https://github.com/JCSDA/saber/tree/develop/test/mains/type_model.F90#L799) | load ensemble data | **model** :  Model - class(model_type) - inout<br>**mpl** :  MPI data - type(mpl_type) - inout<br>**nam** :  Namelist - type(nam_type) - in<br>**filename** :  Filename ('ens1' or 'ens2') - character(len=*) - in |
| subroutine | [model_generate_obs](https://github.com/JCSDA/saber/tree/develop/test/mains/type_model.F90#L864) | generate observations locations | **model** :  Model - class(model_type) - inout<br>**mpl** :  MPI data - type(mpl_type) - inout<br>**nam** :  Namelist - type(nam_type) - in |
