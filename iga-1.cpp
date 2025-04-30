#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>


// B-spline基底関数の計算
double bspline(int i, int p, double u, const std::vector<double>& knots){
    if (p == 0) {
        //0次基底関数(区間内なら1、外なら0)
        if(u>=knots[i]&& u < knots[i+1]){
            return 1.0;
        } else {
            return 0.0;
        }
    }

    //p次基底関数
    
}

double computeBSline(int p, const std::vector<double>& knotVector, const std::vector<double>& controlPoints,int numPoints){

}

int main () {
    //controlPoint
    std::vector<double> controlPoints={0.0, 1.0, 0.0};
    //ノットベクトルの定義
    std::vector<double> knotVector = {0.0, 0.0, 1.0 ,1.0 ,2.0 ,2.0};

    //B-spline曲線を計算
    std::vector<double> curve = computeBSline(1, knotVector,controlPoints,100 );

    std::cout << "B-spline 1D:\n";
    for(int i=0; i< curve.size(); i++){
        std::cout << std::setw(10) << i /100.0 << " : " << curve[i] << std::endl;
    }

    return 0;
}