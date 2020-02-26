/**********
 * Copyright (C) 2020 Xilinx, Inc
 *
 * Licensed under the Apache License, Version 2.0 (the "License"). You may
 * not use this file except in compliance with the License. A copy of the
 * License is located at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
 * WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
 * License for the specific language governing permissions and limitations
 * under the License.
 * 
 **********/

#ifndef __RANSAC_PY_H__
#define __RANSAC_PY_H__

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#ifndef EIGEN_MPL2_ONLY
#define EIGEN_MPL2_ONLY
#endif

namespace fivept_ransac {

namespace py = pybind11;

typedef py::array_t<double, py::array::c_style | py::array::forcecast> pyarray_d;
typedef py::array_t<float, py::array::c_style | py::array::forcecast> pyarray_f;
typedef py::array_t<unsigned int, py::array::c_style | py::array::forcecast> pyarray_ui;

opengv::transformation_t transformationFromArray(pyarray_d &t)
{
    opengv::transformation_t retn = Eigen::Map<Eigen::Matrix<double, 3, 4, Eigen::RowMajor>>(t.mutable_data());;
    return retn;
}

void arrayFromTransformation(pyarray_d &a, opengv::transformation_t &t) 
{
    Eigen::Map<Eigen::Matrix<double, 3, 4, Eigen::RowMajor>>(a.mutable_data()) = t;
}

opengv::essential_t essentialFromArray(double *E)
{
    opengv::essential_t retn = Eigen::Map<Eigen::Matrix<double, 3, 3, Eigen::RowMajor>>(E);
    return retn;
}

} // namespace fivept_ransac

#endif
