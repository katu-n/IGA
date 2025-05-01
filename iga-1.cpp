#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <fstream>

struct Point {
    double x;
    double y;
};

// B-spline基底関数の計算
double bspline(int i, int p, double xi, const std::vector<double>& knots){

    if (p == 0) {
        return (xi >= knots[i] && xi < knots[i + 1]) ? 1.0 : 0.0;
    } else {
        double left = 0.0;
        double right = 0.0;
    
        double denom1 = knots[i + p] - knots[i];
        if (denom1 != 0) {
            left = (xi - knots[i]) / denom1 * bspline(i, p - 1, xi, knots);
        }
        double denom2 = knots[i + p + 1] - knots[i + 1];
        if (denom2 != 0) {
            right = (knots[i + p + 1] - xi) / denom2 * bspline(i + 1, p - 1, xi, knots);
        }
        return left + right;
    } 
}

Point computeBSline(double xi, const std::vector<Point>& controlPoints, const std::vector<double>& knot, int degree) {
    Point result = {0.0, 0.0};

    for(int i = 0; i < controlPoints.size(); ++i) {
        double b = bspline(i, degree, xi, knot);
        result.x += b * controlPoints[i].x;
        result.y += b * controlPoints[i].y;
    }

    return result;

}

int main () {
    //controlPoint
    std::vector<Point> controlPoints={{0,0}, {1,2}, {2,0}, {3,1}, {4,0}};
    //ノットベクトルの定義
    std::vector<double> knot = {0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 4.0};
    int degree = 2; //B-splineの次数

    std::ofstream ofs("bspline.txt");
    //ξはxi
    for(double x1i = knot[degree]; x1i < knot[knot.size()-degree-1]; x1i+=0.01){
        Point pt = computeBSline(x1i, controlPoints, knot, degree);
        ofs << pt.x << " " << pt.y << "\n";
    }

    ofs.close();
    return 0;
}