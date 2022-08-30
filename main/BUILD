load("@rules_cc//cc:defs.bzl", "cc_binary")

cc_binary(
    name = "main",
    srcs = ["main.cc"],
    includes = ["/opt/boost/include"],
    deps = [ 
        "//lib:hh",
        "//lib:neuromorphic",
    ],  
    copts = [
        "-fpic",
    ],
    linkopts = ["-lpthread"],
)
