load("@rules_cc//cc:defs.bzl", "cc_library")

cc_library(
    name = "hh",
    srcs = ["hh.cc"],
    hdrs = ["hh.h"],
    visibility = ["//main:__pkg__"],
)

cc_library(
    name = "neuromorphic",
    srcs = ["neuromorphic.cc"],
    hdrs = ["neuromorphic.h", "hh.h"],
    visibility = ["//main:__pkg__"],
)
