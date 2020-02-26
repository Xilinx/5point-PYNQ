#  Copyright (C) 2020 Xilinx, Inc
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.

from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext as _build_ext
import setuptools
import os
import sys
import pathlib
import pygit2
import distro
from pynq.utils import build_py


__author__ = "Giuseppe Natale"
__copyright__ = "Copyright 2020, Xilinx"
__email__ = "pynq_support@xilinx.com"


# global variables
module_name = "fivepoint_pynq"
data_files = []

# git repos
eigen_git = "https://gitlab.com/libeigen/eigen.git"
opengv_git = "https://github.com/laurentkneip/opengv.git"


class CustomExtension(Extension):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


class _get_pybind_include(object):
    """Helper class to determine the pybind11 include path.

    The purpose of this class is to postpone importing pybind11 until it is
    actually installed, so that the `get_include()` method can be invoked.
    """
    def __init__(self, user=False, pep517=False):
        self.user = user
        self.pep517 = pep517

    def __str__(self):
        # fix to https://github.com/pybind/pybind11/issues/1067
        import pybind11
        if self.pep517:
            # recursing <root>/lib/pythonXX/site-packages/pybind11/__init__.py
            # to get <root>/include
            root_path = pathlib.Path(pybind11.__file__).parent.parent.parent.\
                parent.parent
            include_path = pathlib.Path(root_path, "include")
            for found in include_path.rglob("pybind11.h"):
                return str(found.parent.parent)
        return pybind11.get_include(self.user)


ext_modules = [
    CustomExtension(
        "{}.fivept_ransac".format(module_name),
        [str(pathlib.Path(module_name, "ransac", "ransac_py.cpp"))],
        include_dirs=[
            # Path to pybind11 headers
            _get_pybind_include(),
            _get_pybind_include(user=True),
            _get_pybind_include(pep517=True),
        ],
        libraries=[
            "opengv"
        ],
        language='c++'
    ),
]


# As of Python 3.6, CCompiler has a `has_flag` method.
# cf http://bugs.python.org/issue26689
def _has_flag(compiler, flagname):
    """Return a boolean indicating whether a flag name is supported on
    the specified compiler.
    """
    import tempfile
    with tempfile.NamedTemporaryFile('w', suffix='.cpp') as f:
        f.write('int main (int argc, char **argv) { return 0; }')
        try:
            compiler.compile([f.name], extra_postargs=[flagname])
        except setuptools.distutils.errors.CompileError:
            return False
    return True


def _cpp_flag(compiler):
    """Return the -std=c++[11/14/17] compiler flag.

    The newer version is preferred over c++11 (when it is available).
    """
    flags = ['-std=c++17', '-std=c++14', '-std=c++11']

    for flag in flags:
        if _has_flag(compiler, flag):
            return flag

    raise RuntimeError('Unsupported compiler -- at least C++11 support '
                       'is needed!')


def _get_cmake():
    """Get correct cmake based on OS."""
    dist = distro.id()
    if dist in ["ubuntu", "debian"]:
        return "cmake"
    elif dist in ["rhel", "centos"]:
        return "cmake3"
    else:
        raise OSError("Current OS is not supported.")


class build_ext(_build_ext):
    """A custom build extension for adding compiler-specific options."""
    c_opts = {
        'unix': ["-O3", "-fno-strict-aliasing"],
    }
    l_opts = {
        'unix': ["-O3", "-fno-strict-aliasing"],
    }

    def run(self):
        for ext in self.extensions:
            if type(ext) is CustomExtension:
                # try:
                #     import pyopengv
                # except ImportError:
                #     # pyopengv not found in the system, install locally
                #     self.build_opengv_cmake(ext)
                # To avoid potential issues with linked libraries versions
                # (system libraries vs e.g. libraries installed in the conda
                # environment) we are currently forcing the internal
                # compilation of opengv regardless of whether it is already
                # installed or not
                self.build_opengv_cmake(ext)
        super().run()

    def build_opengv_cmake(self, ext):
        cwd = pathlib.Path().absolute()

        build_temp = pathlib.Path(self.build_temp).absolute()
        pathlib.Path(build_temp, "build").mkdir(parents=True, exist_ok=True)
        extdir = pathlib.Path(self.get_ext_fullpath(ext.name)).absolute()
        pathlib.Path(extdir.parent, "include").mkdir(parents=True,
                                                     exist_ok=True)

        # get Eigen
        eigendir = pathlib.Path(extdir.parent, "include", "eigen")
        pygit2.clone_repository(eigen_git, str(eigendir),
                                checkout_branch="3.3")
        # get OpenGV
        opengv_temp = pathlib.Path(build_temp, "opengv")
        pygit2.clone_repository(opengv_git,
                                str(opengv_temp)).update_submodules(init=True)
        # example of cmake args
        config = 'Debug' if self.debug else 'Release'
        cmake_args = [
            "-DCMAKE_BUILD_TYPE=" + config,
            "-DBUILD_PYTHON=ON",
            "-DEIGEN_INCLUDE_DIR=" + str(eigendir),
            "-DCMAKE_INSTALL_PREFIX=" + str(extdir.parent),
            "-DPYTHON_INSTALL_DIR=" + str(extdir.parent),
            "-DPYTHON_EXECUTABLE=" + sys.executable,
        ]
        build_args = [
            "--config", config,
            "--", "CPPFLAGS=-DEIGEN_MPL2_ONLY"
        ]
        os.chdir(str(pathlib.Path(build_temp, "build")))
        cmake = _get_cmake()
        self.spawn([cmake, str(opengv_temp)] + cmake_args)
        if not self.dry_run:
            self.spawn([cmake, "--build", "."] + build_args)
            self.spawn(["make", "install"])
        os.chdir(str(cwd))
        ext.include_dirs.append(str(eigendir))
        ext.include_dirs.append(str(pathlib.Path(extdir.parent, "include")))
        libdir = str(pathlib.Path(extdir.parent, "lib"))
        ext.library_dirs.append(libdir)
        # statically link libopengv using `extra_objects`
        ext.extra_objects = ['{}/lib{}.a'.format(libdir, "opengv")]

    def build_extensions(self):
        ct = self.compiler.compiler_type
        opts = self.c_opts.get(ct, [])
        link_opts = self.l_opts.get(ct, [])
        if ct == 'unix':
            opts.append('-DVERSION_INFO="{}"'.format(
                self.distribution.get_version()))
            opts.append(_cpp_flag(self.compiler))
            if _has_flag(self.compiler, '-fvisibility=hidden'):
                opts.append('-fvisibility=hidden')
        for ext in self.extensions:
            ext.extra_compile_args = opts
            ext.extra_link_args = link_opts
        super().build_extensions()


def extend_package(path):
    if os.path.isdir(path):
        data_files.extend(
            [os.path.join("..", root, f)
             for root, _, files in os.walk(path) for f in files]
        )
    elif os.path.isfile(path):
        data_files.append(os.path.join("..", path))

with open("README.md", encoding="utf-8") as fh:
    readme_lines = fh.readlines()[4:]
long_description = ("".join(readme_lines))

extend_package(os.path.join(module_name, "notebooks"))
setup(name=module_name,
      version="1.0",
      description="5-point Relative Pose Problem for PYNQ",
      long_description=long_description,
      long_description_content_type="text/markdown",
      author="Giuseppe Natale",
      author_email="pynq_support@xilinx.com",
      url="https://github.com/Xilinx/5point-PYNQ",
      packages=find_packages(),
      download_url="https://github.com/Xilinx/5point-PYNQ",
      ext_modules=ext_modules,
      package_data={
          "": data_files,
      },
      python_requires='>=3.5.2',
      # keeping 'setup_requires' only for readability - relying on
      # pyproject.toml and PEP 517/518
      setup_requires=[
          "pynq>=2.5.1",
          "pybind11>=2.4",
          "pygit2",
          "distro"
      ],
      install_requires=[
          "pynq>=2.5.1",
          "pybind11>=2.4",
          "jupyter",
          "websockets",
          "nest_asyncio"
      ],
      extras_require={
          ':python_version<"3.6"': [
              'ipython==7.9'
          ]
      },
      entry_points={
          "pynq.notebooks": [
              "5point = {}.notebooks".format(module_name)
          ]
      },
      cmdclass={
          "build_py": build_py,
          "build_ext": build_ext,
      },
      zip_safe=False,
      license="Apache License 2.0"
      )
