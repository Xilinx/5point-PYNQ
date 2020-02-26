/*********
 * MIT License
 * 
 * Copyright (c) 2018 NECSTLab, Politecnico di Milano
 * Copyright (c) 2020 Xilinx, Inc
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *********/

#ifndef __RANSAC_H__
#define __RANSAC_H__

#include <opengv/types.hpp>
#include <opengv/triangulation/methods.hpp>
#include <opengv/relative_pose/CentralRelativeAdapter.hpp>
#include <Eigen/StdVector>

#include "ransac_py.h"

#define POINTS_PER_ITER 5
#define EXTRA_RANSAC_POINTS 3
#define TOTAL_POINTS_PER_ITER (POINTS_PER_ITER + EXTRA_RANSAC_POINTS)

namespace fivept_ransac {

namespace py = pybind11;

bool recover_transformation(double *essentials, 
                            unsigned int *essential_sizes,
                            unsigned int *indices,
                            opengv::relative_pose::CentralRelativeAdapter &adapter, 
                            int ransac_iter,
                            opengv::transformation_t &outModel, 
                            opengv::transformation_t &inverseOutModel)
{
    unsigned int tests = essential_sizes[ransac_iter + 1] - essential_sizes[ransac_iter];
    if(tests == 0) {
        return false;
    }
    opengv::transformations_t transformations(tests * 4);
    opengv::transformations_t inverseTransformations(tests * 4);

    //now decompose each essential matrix into transformations and find the
    //right one
    Eigen::Matrix3d W = Eigen::Matrix3d::Zero();
    W(0,1) = -1;
    W(1,0) = 1;
    W(2,2) = 1;

    int indices_offset = ransac_iter * TOTAL_POINTS_PER_ITER;
    unsigned int i;
    size_t o;

    for(i = essential_sizes[ransac_iter], o = 0; 
        i < essential_sizes[ransac_iter + 1]; i++, o++)
    {
        // decompose
        Eigen::MatrixXd tempEssential = essentialFromArray(essentials + i*9);
        Eigen::JacobiSVD< Eigen::MatrixXd > SVD(
                tempEssential,
                Eigen::ComputeFullV | Eigen::ComputeFullU );
        Eigen::VectorXd singularValues = SVD.singularValues();

        // check for bad essential matrix
        if( singularValues[2] > 0.001 ) {};
        // continue; //singularity constraints not applied -> removed because too harsh
        if( singularValues[1] < 0.75 * singularValues[0] ) {};
        // continue; //bad essential matrix -> removed because too harsh

        // maintain scale
        double scale = singularValues[0];

        // get possible rotation and translation vectors
        opengv::rotation_t Ra = SVD.matrixU() * W * SVD.matrixV().transpose();
        opengv::rotation_t Rb = SVD.matrixU() * W.transpose() * SVD.matrixV().transpose();
        opengv::translation_t ta = scale*SVD.matrixU().col(2);
        opengv::translation_t tb = -ta;

        // change sign if det = -1
        if( Ra.determinant() < 0 ) Ra = -Ra;
        if( Rb.determinant() < 0 ) Rb = -Rb;

        //derive transformations
        opengv::transformation_t transformation;
        transformation.col(3) = ta;
        transformation.block<3,3>(0,0) = Ra;
        transformations[o*4 + 0] = transformation;
        transformation.col(3) = ta;
        transformation.block<3,3>(0,0) = Rb;
        transformations[o*4 + 1] = transformation;
        transformation.col(3) = tb;
        transformation.block<3,3>(0,0) = Ra;
        transformations[o*4 + 2] = transformation;
        transformation.col(3) = tb;
        transformation.block<3,3>(0,0) = Rb;
        transformations[o*4 + 3] = transformation;

        // derive inverse transformations
        for(size_t j = 0; j < 4; j++)
        {
            opengv::transformation_t inverseTransformation;
            inverseTransformation.block<3,3>(0,0) =
                    transformations[o*4 + j].block<3,3>(0,0).transpose();
            inverseTransformation.col(3) =
                    -inverseTransformation.block<3,3>(0,0)*transformations[o*4 + j].col(3);
            inverseTransformations[o*4 + j] = inverseTransformation;
        }
    }
    double bestQuality = 1000000.0;
    int bestQualityIndex = -1;

    for(i = 0; i < tests; i++)
    {
        // collect qualities for each of the four solutions solution
        Eigen::Matrix<double,4,1> p_hom;
        p_hom[3] = 1.0;

        for(size_t j = 0; j<4; j++)
        {
            // prepare variables for triangulation and reprojection
            adapter.sett12(transformations[i*4 + j].col(3));
            adapter.setR12(transformations[i*4 + j].block<3,3>(0,0));

            // go through all features and compute quality of reprojection
            double quality = 0.0;

            for( int k = 0; k < TOTAL_POINTS_PER_ITER; k++ )
            {
                p_hom.block<3, 1>(0, 0) =
                        opengv::triangulation::triangulate2(adapter, indices[indices_offset + k]);
                opengv::bearingVector_t reprojection1 = p_hom.block<3, 1>(0, 0);
                opengv::bearingVector_t reprojection2 = inverseTransformations[i*4 + j] * p_hom;
                reprojection1 = reprojection1 / reprojection1.norm();
                reprojection2 = reprojection2 / reprojection2.norm();
                opengv::bearingVector_t f1 = adapter.getBearingVector1(k);
                opengv::bearingVector_t f2 = adapter.getBearingVector2(k);

                // bearing-vector based outlier criterium (select threshold accordingly):
                // 1-(f1'*f2) = 1-cos(alpha) \in [0:2]
                double reprojError1 = 1.0 - (f1.dot(reprojection1));
                double reprojError2 = 1.0 - (f2.dot(reprojection2));

                double error = reprojError1 + reprojError2;
                quality += error;

                if(bestQualityIndex != -1 && quality >= bestQuality) {
                    break;
                }
            }

            // is quality better? (lower)
            if( quality < bestQuality )
            {
                bestQuality = quality;
                bestQualityIndex = i*4 + j;
            }
        }
    }
    if( bestQualityIndex == -1 ) {
        return false; // no solution found
    }

    outModel.col(3) = transformations[bestQualityIndex].col(3);
    outModel.block<3,3>(0,0) = transformations[bestQualityIndex].block<3,3>(0,0);
    inverseOutModel.col(3) = inverseTransformations[bestQualityIndex].col(3);
    inverseOutModel.block<3,3>(0,0) = inverseTransformations[bestQualityIndex].block<3,3>(0,0);
    return true;
}

double compute_quality(opengv::relative_pose::CentralRelativeAdapter &adapter,
                       opengv::transformation_t &model, 
                       opengv::transformation_t &inverseSolution,
                       double ransac_threshold, double best_quality)
{
    opengv::translation_t translation = model.col(3);
    opengv::rotation_t rotation = model.block<3,3>(0,0);
    adapter.sett12(translation);
    adapter.setR12(rotation);

    Eigen::Matrix<double,4,1> p_hom;
    p_hom[3] = 1.0;

    double quality = 0;

    for(size_t i = 0; i < adapter.getNumberCorrespondences(); i++) {
        p_hom.block<3,1>(0,0) =
        opengv::triangulation::triangulate2(adapter, i);
        opengv::bearingVector_t reprojection1 = p_hom.block<3,1>(0,0);
        opengv::bearingVector_t reprojection2 = inverseSolution * p_hom;
        reprojection1 = reprojection1 / reprojection1.norm();
        reprojection2 = reprojection2 / reprojection2.norm();
        opengv::bearingVector_t f1 = adapter.getBearingVector1(i);
        opengv::bearingVector_t f2 = adapter.getBearingVector2(i);

        //bearing-vector based outlier criterium (select threshold accordingly):
        //1-(f1'*f2) = 1-cos(alpha) \in [0:2]
        double reprojError1 = 1.0 - (f1.transpose() * reprojection1);
        double reprojError2 = 1.0 - (f2.transpose() * reprojection2);

        double error = reprojError1 + reprojError2;
        if(error > ransac_threshold) {
            error = 1;
        }

        quality += error;

        if(quality > best_quality) {
            return quality;
        }
    }

    return quality;
}

double ransac(pyarray_d &essentials, pyarray_ui &essential_sizes, int num_ransac_iter,
              pyarray_ui &indices, pyarray_f &p1, pyarray_f &p2, int num_points,
              pyarray_d &transformation, double ransac_threshold)
{
    double best_quality;
    auto E_buf = essentials.request();
    double *_essentials = (double *) E_buf.ptr;
    auto Es_buf = essential_sizes.request();
    unsigned int *_essential_sizes = (unsigned int *) Es_buf.ptr;
    opengv::transformation_t outModel = transformationFromArray(transformation);

    // prepare bearing vectors and central relative adapter
    opengv::bearingVectors_t b1(num_points), b2(num_points);
    opengv::relative_pose::CentralRelativeAdapter adapter(b1, b2);

    auto buf1 = p1.request(), buf2 = p2.request();

    float *ptr1 = (float *) buf1.ptr;
    float *ptr2 = (float *) buf2.ptr;

    for (int j = 0; j < num_points; j++) {
        b1[j] = Eigen::Vector3d(ptr1[j*2], ptr1[j*2 + 1], 1.0);
        b2[j] = Eigen::Vector3d(ptr2[j*2], ptr2[j*2 + 1], 1.0);
        b1[j] = b1[j] / b1[j].norm();
        b2[j] = b2[j] / b2[j].norm();
    }

    auto bufidx = indices.request();
    unsigned int *ptridx = (unsigned int *) bufidx.ptr;

    opengv::transformation_t inverseModel;

    best_quality = num_points * 10;
    int best_iteration = -1;

    for(int i = 0; i < num_ransac_iter; i++) {
        bool result = recover_transformation(_essentials, _essential_sizes,
                                             ptridx, adapter, i, 
                                             outModel, inverseModel);
        if(result) {
            double cur_quality = compute_quality(adapter, outModel, inverseModel, 
                                                 ransac_threshold, best_quality);
            if(cur_quality < best_quality) {
                best_quality = cur_quality;
                best_iteration = i;
            }
        }
    }
    if(best_iteration >= 0) {
        // recompute best transformation
        recover_transformation(_essentials, _essential_sizes, ptridx, adapter, 
                               best_iteration, outModel, inverseModel);
    }
    arrayFromTransformation(transformation, outModel);
    return best_quality;
}

} // namespace fivept_ransac

#endif
