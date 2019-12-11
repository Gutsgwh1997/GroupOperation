#include <Eigen/Core>
#include <Eigen/Geometry>
#include <iostream>
#include "SE3.h"

int main() {
    Eigen::AngleAxisd rotation_vector(M_PI / 4., Vector3d(0, 0, 1));
    Eigen::Quaterniond q(rotation_vector);
    Vector3d t(1, 2, 3);
    Vector3d xyz(1, 2, 1);

    SE3 SE(q, t);
    std::cout << SE << std::endl;

    Vector6d se3 = SE.log();
    std::cout << "The se3 is: " << se3.transpose() << std::endl;

    SE3 SES = SE.exp(se3);
    std::cout << "Then exp() is: \n" << SES << std::endl;

    Vector3d after = SE.map(xyz);
    std::cout << "Point before transform: " << xyz.transpose() << std::endl;
    std::cout << "Point after transform: " << after.transpose() << std::endl;

    SE3 SE_SES = SE*SES;
    std::cout<<"Valide operator*: \n"<<SE_SES<<std::endl;


    return 0;
}
