#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <fstream>

struct coordinate
{
    double x,y;
};

class Bspline
{
    public:
        std::vector<double> knot;
        std::vector<coordinate> CP;
        int p; // degree
        double newKnot;
        int step= 100;
        Bspline(int p,const std::vector<coordinate>& CP):p(p),CP(CP)
        {
            int ncp = CP.size();
            int m   = ncp + p + 1;
            knot.resize(m);

            for(int i=0; i<=p; i++) 
           {
                knot[i] = 0.0;
                knot[m-1-i] = 1.0;
           } 

           int interKnot = m -2*(p+1);
           if(interKnot > 0)
           {
                double stepKnot = 1.0/(interKnot+1);
                for(int i=0; i<interKnot; i++)
                {
                    knot[p+1+i] = stepKnot*(i+1);
                }
           }
        }
        ~Bspline(){}

        double basisFunction(int i, int degree, double t);
        void KnotInsertion(double newKnot);
        void calCoordinate();
        void outputBasisFunction();
        void printCSV(std::ofstream& ofs,double xi,double x, double y);
};

double Bspline::basisFunction(int i,int degree, double t)
{
    if(degree==0)
    {
        if((knot[i]<=t && t<knot[i+1]) || (std::abs(t-1.0)<1e-9 && i == knot.size()- this->p -2))
        {
            return 1.0;
        }
        return 0.0;
    }

    double term1 =0.0, term2=0.0;
    double denom1 = knot[i+degree]-knot[i];
    double denom2 = knot[i+degree+1]-knot[i+1];

    if(denom1 > 1e-9)
    {
        term1 = (t - knot[i])/denom1 * basisFunction(i, degree-1, t);
    }

    if(denom2 > 1e-9)
    {
        term2 = (knot[i+degree+1]-t)/denom2 * basisFunction(i+1, degree-1, t);
    }

    return term1+term2;
}

void Bspline::printCSV(std::ofstream& ofs, double xi, double x, double y)
{
        ofs << std::fixed << std::setprecision(4) << xi << ","
            << std::setprecision(6) << x << "," << y << std::endl;
}


void Bspline::calCoordinate()
{
    std::ofstream ofs("b_splineCurve_KnotInsert.csv");

    ofs << "t,x,y" << std::endl;

        for(int i=0; i<=step; i++)
        {
            double xi = static_cast<double>(i)/step;

            double x=0.0, y=0.0;
            for(int j=0; j<CP.size(); j++)
            {
                double N = basisFunction(j,p,xi);
                x += CP[j].x*N;
                y += CP[j].y*N;
            }
            printCSV(ofs,xi,x,y);
        }
}

void Bspline::KnotInsertion(double newKnot)
{
    int k = 0;

    for(size_t i=0; i<knot.size()-1; i++)
    {
        if(knot[i] <= newKnot && newKnot < knot[i+1])
        {
            k = i;
            break;
        }
    }


        printf("k=%d \n",k);

    std::vector<coordinate> newCP;
    int n = CP.size();

    newCP.resize(n+1);

    for(int i=0; i<=n; i++)
    {
        double alpha;
        if(i <= k-p)
        {
            alpha = 1;
        }
        else if(i>=k-p+1 && i<=k)
        {
            alpha = (newKnot - knot[i])/(knot[i+p]-knot[i]);
        }
        else
        {
            alpha = 0;
        }
        coordinate tmp;

        newCP[i].x = alpha*CP[i].x + (1.0-alpha)*CP[i-1].x;
        newCP[i].y = alpha*CP[i].y + (1.0-alpha)*CP[i-1].y;

        printf("alpha=%f tmp.y=%f \n", alpha,tmp.y);
    }
    CP = newCP;
    knot.insert(knot.begin()+k+1,newKnot);
}

void Bspline::outputBasisFunction()
{
    std::ofstream ofsB("knotIncert_basisFunction.csv");

    for(size_t j=0; j<CP.size(); j++)
    {
        ofsB << "N" << j << ",";
    }
    ofsB <<  std::endl;

        for(int i=0; i<=step; i++)
        {
            double xi = static_cast<double>(i)/step;

            for(int j=0; j<CP.size(); j++)
            {
                double N = basisFunction(j,p,xi);
                ofsB <<std::fixed << std::setprecision(6) << N << ",";
            }
            ofsB << std::endl;
        }
}

int main()
{
    int p = 3; // degree
    std::vector<coordinate> cps = {
        {0.0, 0.0},
        {0.0, 2.0},
        {1.0, 2.0},
        {2.0, 0.5},
        {5.0, -4.0},
        {5.0, 2.0},
        {12.0, 0.5}
    };

    double newKnot=0.9;

    Bspline bs(p,cps);
    bs.KnotInsertion(newKnot);
    bs.calCoordinate();
    bs.outputBasisFunction();

    std::cout << "Knot Vector:";
    for(const auto& k : bs.knot)
    {
        std::cout << std::fixed << std::setprecision(4) << k << " ";
    }
    std::cout << std::endl;

    return 0;
}
