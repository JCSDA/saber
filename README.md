### Continuous integration:
| Platform      |  JCSDA-internal       | JCSDA      |
| ------------- | ------------- |------------- |
| GNU           | [![AWS-gnu](https://codebuild.us-east-1.amazonaws.com/badges?uuid=eyJlbmNyeXB0ZWREYXRhIjoieXA5WFFUZk5NMDVvV0ZkZVBQUDRVeFN3VDk5aVkvZHJ0K3ZWaEl6RlVZNGdCTEI0Y283QU5TTzVTS0k0N0hNYjl4allwY21SRlFWVjJYTEFjSlJUUlZVPSIsIml2UGFyYW1ldGVyU3BlYyI6ImxZNTZLc3VXcGNVYktCeVQiLCJtYXRlcmlhbFNldFNlcmlhbCI6MX0%3D&branch=develop)](https://console.aws.amazon.com/codesuite/codebuild/469205354006/projects/saber-internal-gnu/history) | [![AWS-gnu](https://codebuild.us-east-1.amazonaws.com/badges?uuid=eyJlbmNyeXB0ZWREYXRhIjoiV2dVMmxFVENKL2dCVzN5UlgyZHJuSmhvbTV6dDhOalYwTEJDaXdZWGFDbXp2YlU4VzdsV3ZRNm9mT25mRnM3NlVYWXE2R2pmYVlZbWhxbHJ1OXFpdzVjPSIsIml2UGFyYW1ldGVyU3BlYyI6Ilp2T04vNnBRR0xFYmQ3UzAiLCJtYXRlcmlhbFNldFNlcmlhbCI6MX0%3D&branch=develop)](https://console.aws.amazon.com/codesuite/codebuild/469205354006/projects/automated-testing-saber-gnu/history)
| Intel         | [![AWS-intel](https://codebuild.us-east-1.amazonaws.com/badges?uuid=eyJlbmNyeXB0ZWREYXRhIjoidC9ZWmlyNU8xZmdjd0kxbFJzcGVNTlhJSDdBcFJ4RUdwNjNmcnFzQ1VWUUNaMWFEZkwvbHlkZUxTaTZIZlQyWWxOMGtvVzRaTlpRNGdjbFVUK0ZaRDFvPSIsIml2UGFyYW1ldGVyU3BlYyI6IllwQlZTb2JNdnJjOEo5TlgiLCJtYXRlcmlhbFNldFNlcmlhbCI6MX0%3D&branch=develop)](https://console.aws.amazon.com/codesuite/codebuild/469205354006/projects/saber-internal-intel/history) | [![AWS-intel](https://codebuild.us-east-1.amazonaws.com/badges?uuid=eyJlbmNyeXB0ZWREYXRhIjoiYUROTE5DZVdranpBQTBKbTlBam1vb2pVWXJteDdEMk1RLzhWdmlQU2NUQUhueFF2UnhINWxDcGZ1eWFqcFpBUVRDMGpYdVhzSWdmazNYcmRDeUdOd0xRPSIsIml2UGFyYW1ldGVyU3BlYyI6IjhqZnUxOHpObWFGSnFtUzYiLCJtYXRlcmlhbFNldFNlcmlhbCI6MX0%3D&branch=develop)](https://console.aws.amazon.com/codesuite/codebuild/469205354006/projects/automated-testing-saber-intel/history)
| CLANG         | [![AWS-clang](https://codebuild.us-east-1.amazonaws.com/badges?uuid=eyJlbmNyeXB0ZWREYXRhIjoicnRqWEh6YUg1UEx2OWE5QVpXb2RjVDFCeitJV3ROaEkxVGVnYnRNYWMzR0J0Z2xPZFhTZlEvVUFiL1BoUjJzcVh3V3BSaTRaSVFnK2dSdGtMcnd5S2o4PSIsIml2UGFyYW1ldGVyU3BlYyI6IjFVTEtZRTNpQXJMR0NYRCsiLCJtYXRlcmlhbFNldFNlcmlhbCI6MX0%3D&branch=develop)](https://console.aws.amazon.com/codesuite/codebuild/469205354006/projects/saber-internal-clang/history) | [![AWS-clang](https://codebuild.us-east-1.amazonaws.com/badges?uuid=eyJlbmNyeXB0ZWREYXRhIjoiL3NrZ05zdXQzbmlhOTJOT0RVanBwKzhocXhIb0tpdnFFMzAzdjd6RmN4V0FpRTJMVkdYcGJoVS9CTlE0L3dXS3JvclZxZU12U0lVWjdBb3krZ2xzODBBPSIsIml2UGFyYW1ldGVyU3BlYyI6IklHcGQ0VUJNOWdzNHNyWE0iLCJtYXRlcmlhbFNldFNlcmlhbCI6MX0%3D&branch=develop)](https://console.aws.amazon.com/codesuite/codebuild/469205354006/projects/automated-testing-saber-clang/history)
| Code Coverage | [![codecov](https://codecov.io/gh/JCSDA-internal/saber/branch/develop/graph/badge.svg?token=GKZ5TMF2GW)](https://codecov.io/gh/JCSDA-internal/saber) |

# SABER
&copy; Copyright 2019 UCAR

This software is licensed under the terms of the Apache Licence Version 2.0
which can be obtained at [http://www.apache.org/licenses/LICENSE-2.0](http://www.apache.org/licenses/LICENSE-2.0).

The compilation of SABER is made easier by using the [saber-bundle](CI/README.md).

Documentation can be found on the [JEDI Documentation website](https://jointcenterforsatellitedataassimilation-jedi-docs.readthedocs-hosted.com/en/latest/index.html).

## BUMP

The BUMP (Background error on an Unstructured Mesh Package) library estimates and applies background error covariance-related operators, defined on an unstructured mesh.

Most of the BUMP code is distributed under the [CeCILL-C license](http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html) (Copyright 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT).

The fact that you are downloading this code means that you have had knowledge of the [CeCILL-C license](http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html) and that you accept its terms.

## Other libraries coming soon ...
