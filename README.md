# Testing Data Set

Our test suite drew from three distinct benchmarks: **NPB v3.3**, **PB Suite V4.2**, and the **PCB Benchmark**. For our experiments, we measured the performance of the applications using input **Class B** for the NPB and the **LARGE_DATASET** for the PB and PCB suites.

![Test Suite Image](image.png/image.png)

Recall that we are focusing on challenge patterns for optimizing compilers. We have selected specific subroutines from the NPB and PB suites that contain such cases, resulting in a total of **16 experimental sections**. For the PCB, we developed **28 use cases**, each representing one of the six challenge scenarios. In total, our testing dataset contains **44 experimental sections**. 

**Table 2** below displays our test suite.

| Benchmark | Suite Version | Input Dataset    | Experimental Sections |
|-----------|---------------|-------------------|-----------------------|
| NPB       | v3.3          | Class B           | 12                    |
| PB        | V4.2          | LARGE_DATASET     | 4                    |
| PCB       | -             | LARGE_DATASET     | 28                    |
| **Total** |               |                   | **44**                |

*Table 2: Overview of the Test Suite*

# TESTING_SUITE_AOTs
