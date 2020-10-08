# Module tools_samp

| Type | Name | Purpose | Arguments |     | Type | Intent |
| :--: | :--: | :------ | ----: | :-------- | :--: | :----: |
| subroutine | [initialize_sampling](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/tools_samp.F90#L54) | intialize sampling | **mpl**<br>**rng**<br>**area**<br>**n_loc**<br>**lon_loc(n_loc)**<br>**lat_loc(n_loc)**<br>**mask_loc(n_loc)**<br>**rh_loc(n_loc)**<br>**loc_to_glb(n_loc)**<br>**ntry**<br>**nrep**<br>**ns2_glb**<br>**sam2_glb(ns2_glb)**<br>**fast**<br>**verbosity**<br>**n_uni**<br>**uni_to_loc(:)**<br>**tree_uni** |  MPI data<br> Random number generator<br> Global domain area<br> Number of points (local)<br> Longitudes (local)<br> Latitudes (local)<br> Mask (local)<br> Horizontal support radius (local)<br> Local to global index<br> Number of tries<br> Number of replacements<br> Number of samplings points (global)<br> Horizontal sampling index (global)<br> Fast sampling flag<br> Verbosity flag<br> Universe size<br> Universe to local index<br> Universe KD-tree | type(mpl_type)<br>type(rng_type)<br>real(kind_real)<br>integer<br>real(kind_real)<br>real(kind_real)<br>logical<br>real(kind_real)<br>integer<br>integer<br>integer<br>integer<br>integer<br>logical<br>logical<br>integer<br>integer<br>type(tree_type) | inout<br>inout<br>in<br>in<br>in<br>in<br>in<br>in<br>in<br>in<br>in<br>in<br>out<br>in<br>in<br>in<br>in<br>in |
| subroutine | [initialize_sampling_local](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/tools_samp.F90#L192) | intialize sampling, local, based on ATLAS octahedral grid | **mpl**<br>**area**<br>**n_loc**<br>**n_loc_eff**<br>**n_glb_eff**<br>**lon_loc(n_loc)**<br>**lat_loc(n_loc)**<br>**mask_loc(n_loc)**<br>**rh_loc(n_loc)**<br>**loc_to_glb(n_loc)**<br>**ns2_glb**<br>**ns1_glb_eff**<br>**lon1_glb_eff(:)**<br>**lat1_glb_eff(:)**<br>**rh1_glb_eff(:)**<br>**sam1_glb_eff(:)**<br>**verbosity**<br>**n_uni**<br>**uni_to_loc(:)**<br>**tree_uni** |  MPI data<br> Global domain area<br> Number of points (local)<br> Number of points (local, effective)<br> Number of points (global, effective)<br> Longitudes (local)<br> Latitudes (local)<br> Mask (local)<br> Horizontal support radius (local)<br> Local to global index<br> Number of samplings points (global)<br>integer,intent(out) :: ns1_glb_eff<br>real(kind_real),allocatable,intent(out) :: lon1_glb_eff(:)<br>real(kind_real),allocatable,intent(out) :: lat1_glb_eff(:)<br>real(kind_real),allocatable,intent(out) :: rh1_glb_eff(:)<br>integer,allocatable,intent(out) :: sam1_glb_eff(:)<br> Verbosity flag<br> Universe size<br> Universe to local index<br> Universe KD-tree | type(mpl_type)<br>real(kind_real)<br>integer<br>integer<br>integer<br>real(kind_real)<br>real(kind_real)<br>logical<br>real(kind_real)<br>integer<br>integer<br>integer<br>real(kind_real)<br>real(kind_real)<br>real(kind_real)<br>integer<br>logical<br>integer<br>integer<br>type(tree_type) | inout<br>in<br>in<br>in<br>in<br>in<br>in<br>in<br>in<br>in<br>in<br>out<br>out<br>out<br>out<br>out<br>in<br>in<br>in<br>in |
| subroutine | [initialize_sampling_global](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/tools_samp.F90#L528) | intialize sampling, global | **mpl**<br>**rng**<br>**ns1_glb_eff**<br>**lon1_glb_eff(ns1_glb_eff)**<br>**lat1_glb_eff(ns1_glb_eff)**<br>**rh1_glb_eff(ns1_glb_eff)**<br>**sam1_glb_eff(ns1_glb_eff)**<br>**ntry**<br>**nrep**<br>**ns2_glb**<br>**sam2_glb(ns2_glb)**<br>**fast**<br>**verbosity** |  MPI data<br> Random number generator<br>integer,intent(in) :: ns1_glb_eff<br>real(kind_real),intent(in) :: lon1_glb_eff(ns1_glb_eff)<br>real(kind_real),intent(in) :: lat1_glb_eff(ns1_glb_eff)<br>real(kind_real),intent(in) :: rh1_glb_eff(ns1_glb_eff)<br>integer,intent(in) :: sam1_glb_eff(ns1_glb_eff)<br> Number of tries<br> Number of replacements<br> Number of samplings points (global)<br> Horizontal sampling index (global)<br> Fast sampling flag<br> Verbosity flag | type(mpl_type)<br>type(rng_type)<br>integer<br>real(kind_real)<br>real(kind_real)<br>real(kind_real)<br>integer<br>integer<br>integer<br>integer<br>integer<br>logical<br>logical | inout<br>inout<br>in<br>in<br>in<br>in<br>in<br>in<br>in<br>in<br>out<br>in<br>in |