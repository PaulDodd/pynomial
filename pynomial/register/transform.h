#pragma once
// stdlib include
#include <iostream>
#include <vector>
#include <random>
// eigen include
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/Geometry>


namespace pynomial{
    namespace _register{

        template<class Precision>
        Eigen::Quaternion<Precision> QuatFromRotationMatrix(const Eigen::Matrix<Precision,3,3> mat)
        {
            Eigen::Quaternion<Precision> q(mat);
            q.normalize();
            return q;
        }

        Eigen::VectorXd make_quaternion_from_rotation_matrix(const Eigen::MatrixXd& Rotation)
        {
            Eigen::VectorXd output = Eigen::VectorXd::Zero(4);
            const Eigen::MatrixXd& input = Rotation;
            Eigen::Matrix<double,3,3> rot = Eigen::Matrix<double,3,3>::Zero();

            if(input.rows() == input.cols() && (input.rows() == 2 || input.rows()==3) ) // must be square and in 2 or 3 dimensions
            {
                if(input.rows() == 2) // convert 2D rotation to 3D rotaion about z-axis
                {
                    for(int r = 0; r < input.rows(); r++)
                        for(int c = 0; c < input.cols(); c++)
                            rot(r,c) = input(r,c);
                    rot(2,2) = 1.0;
                }
                else
                {
                    rot = input;
                }
                Eigen::Quaternion<double> quat = QuatFromRotationMatrix(rot);
                output[0] = quat.w();
                output[1] = quat.x();
                output[2] = quat.y();
                output[3] = quat.z();
            }
            //TODO: print some error message otherwise.
            return output;
        }
    }
}
