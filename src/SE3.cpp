#include "SE3.h"

/**
 * @brief 一个旋转可以用两个互为相反数的单位四元数表示
 */
void SE3::normalizeRotation(){
    if (_q.w()<0){
        _q.coeffs()*=-1;
    }

    _q.normalize();
}

/**
 * @brief 将三维坐标点xyz按照当前的_q,_t变换
 *
 * @param xyz 三维坐标点
 *
 * @return  变换后的三维坐标点
 */
Vector3d SE3::map(const Vector3d& xyz) const {
    return _q*xyz + _t;
}

/**
 * @brief 李群的对数映射log(),SLAM书上P71页
 *
 * @return 前三维是平移，后三维是旋转
 */
Vector6d SE3::log() const {
    Vector6d res;
    Eigen::Matrix3d _R = _q.toRotationMatrix();
    double d = 0.5 * (_R(0, 0) + _R(1, 1) + _R(2, 2) - 1);   //对应求解theta = arccos(d)
    Vector3d omega;                                          //李群对应的李代数，6维，此变量代表旋转
    Vector3d upsilon;                                        //李群对应的李代数，6维，此变量代表平移

    Vector3d dR = deltaR(_R);
    Eigen::Matrix3d J_inv;

    // 当d接近1时，theta接近0度，自己求解极限，增加数值计算的稳定性
    if (d > 0.99999) {
        omega = 0.5 * dR;
        Eigen::Matrix3d Omega = skew(omega);
        J_inv = Eigen::Matrix3d::Identity() - 0.5 * Omega + (1. / 12.) * (Omega * Omega);
    }
    else{
        double theta = std::acos(d);
        omega = theta / (2 * std::sqrt(1 - d * d)) * dR;  //对应的旋转向量,大小*方向
        Eigen::Matrix3d Omega = skew(omega);
        J_inv = Eigen::Matrix3d::Identity() - 0.5 * Omega //公式(4.26)的逆
            + (1 - theta / (2 * std::tan(theta / 2))) / (theta * theta) * (Omega * Omega);
    }

    upsilon = J_inv*_t;
    for (int i = 0; i < 3; ++i) {
        res[i] = upsilon[i];
    }
    for (int i = 0; i < 3; ++i) {
        res[i + 3] = omega[i];
    }
    return res;
}

/**
 * @brief 李群的指数映射，是static成员函数
 *
 * @param v 输入的6维度李代数，前3维是平移，后三维是旋转
 *
 * @return 返回李群类
 */
SE3 SE3::exp(const Vector6d& v) {
    Vector3d omega;
    for (int i = 0; i < 3; i++)
        omega[i] = v[i + 3];
    Vector3d upsilon;
    for (int i = 0; i < 3; ++i)
        upsilon[i] = v[i];

    double theta = omega.norm();
    Eigen::Matrix3d Omega = skew(omega);

    Eigen::Matrix3d R;
    Eigen::Matrix3d J;
    if(theta < 0.00001){
        R = Eigen::Matrix3d::Identity() + Omega + 0.5 * Omega * Omega;
        J = Eigen::Matrix3d::Identity() + 0.5 * Omega + (1. / 3.) * Omega * Omega;
    }
    else{
        R = Eigen::Matrix3d::Identity() + std::sin(theta) / theta * Omega 
            + (1 - std::cos(theta)) / (theta * theta) * (Omega * Omega);
        J = Eigen::Matrix3d::Identity() + (1 - std::cos(theta)) / (std::pow(theta, 2)) * Omega 
            + (1 - std::sin(theta) / theta) / (theta * theta) * (Omega * Omega);
    }
    return SE3(Eigen::Quaterniond(R), J*upsilon);
}

/**
 * @brief 重载为李群上的乘法.
 *
 * @param tr2 SE3_1*SE3_2中的SE3_2
 *
 * @return 乘法的结果
 */
SE3 SE3::operator*(const SE3& tr2) const{
    SE3 result(*this);   //使用默认复制构造函数，不会有问题
    result._t += _q*tr2._t ;
    result._q *= tr2._q;
    result.normalizeRotation();
    return result;
} 

std::ostream& operator<<(std::ostream& os, const SE3& SE){
    os<<"Rotation is (quaternion w x y z): "<<SE._q.coeffs().transpose()<<std::endl;
    os<<"Translation is: "<<SE._t.transpose();
    return os;
}
