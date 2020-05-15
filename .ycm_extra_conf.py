import os
import sysconfig

import ycm_core

flags = [
    "-x",
    "c++",
    f"-I{os.path.dirname(os.path.abspath(__file__))}",
    "-I/usr/local/include",
    "-Wall",
    "-Werr",
    "-Wconversion",
    "-Weffc++",
    "-std=c++14",
]


def Settings(**kwargs):
    if kwargs["language"] == "cfamily":
        return {"flags": flags}
    if kwargs["language"] == "python":
        pass
