from setuptools import setup, Extension
import pybind11

module = Extension(
    "qpproblem_smo_solver_module",
    sources=["qpproblem_smo_solver_module/src/main.cpp"],
    include_dirs=[pybind11.get_include(), "qpproblem_smo_solver_module/src/"],
    language="c++",
)

setup(
    name="qpproblem_smo_solver",
    version="0.1",
    description="C++ module solving QP problem appearing while training Support Vector Machine (binary classification) with regularization parametr C",
    ext_modules=[module],
)
