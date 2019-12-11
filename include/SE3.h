#ifndef MY_SE3_H
#define MY_SE3_H

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <cmath>
#include <iostream>

// 存储顺序是按照列存储
typedef Eigen::Matrix<double, 6, 1, Eigen::ColMajor> Vector6d;
typedef Eigen::Matrix<double, 3, 1, Eigen::ColMajor> Vector3d;

/**
 * @brief 向量转换为反对称矩阵
 *
 * @param v 三维输入向量
 *
 * @return  三维输出矩阵
 */
Eigen::Matrix3d skew( const Vector3d& v){
    Eigen::Matrix3d m;
    m.fill(0.);
    m(0, 1) = -v(2);
    m(0, 2) = v(1);
    m(1, 2) = -v(0);
    m(1, 0) = v(2);
    m(2, 0) = -v(1);
    m(2, 1) = v(0);
    return m;
}

/**
 * @brief 旋转矩阵转换为旋转向量，求解旋转轴时n = 1/(2*sin(theta))*delta(R),P49
 *
 * @param R 旋转矩阵
 *
 * @return [r32-r23,r13-r31,r21-r12]'
 */
Vector3d deltaR(const Eigen::Matrix3d& R) {
    Vector3d v;
    v(0) = R(2, 1) - R(1, 2);
    v(1) = R(0, 2) - R(2, 0);
    v(2) = R(1, 0) - R(0, 1);
    return v;
}


/**
 * @brief 特殊欧式群
 */
class SE3 {
   public:
       //使用单位四元数与平移向量表示李代数
        SE3(Eigen::Quaterniond q, Vector3d t) : _q(q), _t(t) { normalizeRotation(); }
        //对四元数单位化
        void normalizeRotation();
        // 对xyz做变换
        Vector3d map(const Vector3d& xyz) const;
        // 李群的对数映射
        Vector6d log() const;
        // 指数映射
        static SE3 exp(const Vector6d& v);
        // 李群乘法
        SE3 operator* (const SE3& tr2) const;
        // 显示
        friend std::ostream& operator<<(std::ostream& os, const SE3& SE);

   private:
        Eigen::Quaterniond _q;
        Vector3d _t;
};

#endif
