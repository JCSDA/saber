[![travis develop](https://travis-ci.com/JCSDA/saber.svg?token=zswWHqwVimHTBAygfenZ&branch=develop&logo=travis)](https://travis-ci.com/JCSDA/saber)
[![codecov](https://codecov.io/gh/JCSDA/saber/branch/develop/graph/badge.svg?token=aLmdMnzx1C)](https://codecov.io/gh/JCSDA/saber)
[![AWS-gnu](https://codebuild.us-east-1.amazonaws.com/badges?uuid=eyJlbmNyeXB0ZWREYXRhIjoiV2dVMmxFVENKL2dCVzN5UlgyZHJuSmhvbTV6dDhOalYwTEJDaXdZWGFDbXp2YlU4VzdsV3ZRNm9mT25mRnM3NlVYWXE2R2pmYVlZbWhxbHJ1OXFpdzVjPSIsIml2UGFyYW1ldGVyU3BlYyI6Ilp2T04vNnBRR0xFYmQ3UzAiLCJtYXRlcmlhbFNldFNlcmlhbCI6MX0%3D&branch=develop)](https://us-east-1.console.aws.amazon.com/codesuite/codebuild/projects/automated-testing-saber-gnu/history)
[![AWS-intel](https://codebuild.us-east-1.amazonaws.com/badges?uuid=eyJlbmNyeXB0ZWREYXRhIjoiYUROTE5DZVdranpBQTBKbTlBam1vb2pVWXJteDdEMk1RLzhWdmlQU2NUQUhueFF2UnhINWxDcGZ1eWFqcFpBUVRDMGpYdVhzSWdmazNYcmRDeUdOd0xRPSIsIml2UGFyYW1ldGVyU3BlYyI6IjhqZnUxOHpObWFGSnFtUzYiLCJtYXRlcmlhbFNldFNlcmlhbCI6MX0%3D&branch=develop)](https://us-east-1.console.aws.amazon.com/codesuite/codebuild/projects/automated-testing-saber-intel/history?region=us-east-1)

# SABER
&copy; Copyright 2019 UCAR

This software is licensed under the terms of the Apache Licence Version 2.0
which can be obtained at [http://www.apache.org/licenses/LICENSE-2.0](http://www.apache.org/licenses/LICENSE-2.0).

The compilation of SABER is made easier by using the [saber-bundle](https://github.com/JCSDA/saber-bundle).

## BUMP

The BUMP (B matrix on an Unstructured Mesh Package) library estimates and applies background error covariance-related operators, defined on an unstructured mesh.

Most of the BUMP code is distributed under the [CeCILL-C license](http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html) (Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT).

The fact that you are downloading this code means that you have had knowledge of the [CeCILL-C license](http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html) and that you accept its terms.

Theoretical documentation:
 - about covariance filtering: [covariance_filtering.pdf](doc/bump/pdf/covariance_filtering.pdf)
 - about the NICAS method: [nicas.pdf](doc/bump/pdf/nicas.pdf)
 - about multivariate localization: [multivariate_localization.pdf](doc/bump/pdf/multivariate_localization.pdf)
 - about diffusion and the Matern function: [diffusion_matern_function.pdf](doc/bump/pdf/diffusion_matern_function.pdf)

Code documentation:
 - [Code size and characteristics](doc/bump/CLOC_REPORT.md)
 - [Standalone or online usage](doc/bump/standalone_or_online_usage.md)
 - [Code architecture](doc/bump/code_architecture.md)
 - [Code auto-documentation](doc/bump/code_autodoc.md)
 - [Input data](doc/bump/input_data.md)
 - [Running the code](doc/bump/running_the_code.md)
 - [NCL plots](doc/bump/ncl_plots.md)
 - [Test](doc/bump/test.md)
 - [Adding a new model](doc/bump/adding_a_new_model.md)

## Other libraries coming soon ...
