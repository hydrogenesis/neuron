#!/bin/bash

#bazel build //main:main && bazel-bin/main/main && python plot.py
#bazel build //main:main && bazel-bin/main/main 1 data/ADM21-7sec7_a.CNG.swc && python plot.py
#bazel build //main:main && bazel-bin/main/main 1 data/ADM21-7sec10.CNG.swc && python plot.py
#bazel build //main:main && bazel-bin/main/main 1 data/10_2REDO-850-GM18-Ctl-Ctl-Chow-BNL16A-CA1_Finished2d.CNG.swc && python plot.py
#bazel build //main:main && bazel-bin/main/main 1 data/10_2REDO-850-GM18-Ctl-Ctl-Chow-BNL16A-CA1_Finished2h.CNG.swc && python plot.py
#bazel build //main:main --verbose_failures && bazel-bin/main/main 94 data/2013_03_06_cell03_789_H41_03.CNG.swc && python plot.py
bazel build //main:main --verbose_failures && bazel-bin/main/main 94 data/Ball_and_stick.swc && python plot.py
#bazel build //main:main --compilation_mode=dbg -s && gdb --args bazel-bin/main/main 40 data/ADM21-7sec7_a.CNG.swc
