#!/bin/bash

bazel build //main:main && bazel-bin/main/main && python plot.py
