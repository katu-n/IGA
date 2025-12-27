#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <fstream>

class Bspline
{
    double *B;
    double bnom,term1,term2,coefficient;
    public:
    int p=10;
    double step= 100;
    Bspline()
    {
        B = new double[p+1];
    }
    ~Bspline()
    {
        delete[] B;
    }
    double combination(int a, int p);
    void Bernstein(double x);
    void print(double x);
    void printCSV(std::ofstream& ofs,double x);
};

double Bspline::combination(int a, int p)
{
    if(a==0 || a==p) return 1.0;
    if(a<1 || p < a) return 0;
    if(a > p/2) a = p-a;

    coefficient = 1.0;
    for(int i=1; i<=a; i++)
    {
        coefficient = coefficient*(p-i+1)/i;
    }
    return coefficient;
}

void Bspline::Bernstein(double x)
{
    for(int a=0; a<=p; a++)
    {
        bnom = combination(a,p);

        term1 = std::pow(x, a);
        term2 = std::pow(1.0-x, p-a);

        B[a] = bnom*term1*term2;
    }
}

// void Bspline::print(int x)
// {
//     double sum = 0.0;
//     std::cout << "B-spline Basis Functions (p=" << p << ", x=" << x << "):\n";
//     for(int i=0; i<=p; i++)
//     {
//         sum += B[i];
//         std::cout << "B[" << i << "] = "  << B[i] << std::endl;
//     }
//     std::cout << "SUM = " << sum << std::endl;
// }

void Bspline::printCSV(std::ofstream& ofs, double x)
{
        ofs << std::fixed << std::setprecision(4) << x;
        for(int i=0; i<=p; i++)
        {
            ofs << "," << std::fixed << std::setprecision(6)<< B[i];
        }
        ofs << std::endl;
}

int main()
{
    Bspline bs;

    std::ofstream ofs("bernstein.csv");

    if(!ofs)
    {
        std::cerr << "Error opening csv" << std::endl;
        return -1;
    }

    ofs << "x";
    for(int i=0; i<= bs.p; i++) 
    {
        ofs << ",B" << i;
    }
    ofs << std::endl;

    for(int i=0; i<=bs.step; i++)
    {
        double x = static_cast<double>(i)/bs.step;
        bs.Bernstein(x);
        // bs.print(x);
        bs.printCSV(ofs,x);
    }
    return 0;
}
