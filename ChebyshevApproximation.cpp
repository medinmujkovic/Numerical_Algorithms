#include <iostream>
#include <vector>
#include <cmath>

double pi=4*atan(1);

class ChebyshevApproximation
{
private:
    double xmin,xmax;
    int n,m;
    std::vector<double>c;
    ChebyshevApproximation(std::vector<double>c, double xmin, double xmax): c(c),xmin(xmin),xmax(xmax),m(c.size()-1){}

public:
template <typename FunType>
    ChebyshevApproximation(FunType f, double xmin, double xmax, int n);
    void set_m(int m);
    void trunc(double eps);
    double operator()(double x) const;
    double derivative(double x) const;
    ChebyshevApproximation derivative() const;
    ChebyshevApproximation antiderivative() const;
    double integrate(double a, double b) const;
    double integrate() const;
};

template <typename FunType>
ChebyshevApproximation::ChebyshevApproximation(FunType f, double xmin, double xmax, int n)
    :xmin(xmin),xmax(xmax),n(n){
    m=n; std::vector<double>w,v;
    c.resize(n+1);

    if(xmin>=xmax || n<1) throw std::domain_error("Bad parameters");

    for(int i=0;i<=4*n+3;i++)
        w.push_back(std::cos(pi*i/(2*n+2)));

    for(int i=0;i<=n;i++)
        v.push_back(f((xmin+xmax+(xmax-xmin)*w[2*i+1])/2));
        
    for(int k=0;k<=n;k++)
    {
        double s=0;
        for(int i=0;i<=n;i++)
            s+=v[i]*w[(k*(2*i+1))%(4*n+4)];
        c[k]=2*s/(n+1);
    }
}

void ChebyshevApproximation::set_m(int m)
{
    if(m<=1||m>=n) throw std::domain_error("Bad order");
    this->m=m;
}

void ChebyshevApproximation::trunc(double eps)
{
    if(eps<0) throw std::domain_error("Bad tolerance");
    for(int k=m;k>=0;k--)
    {
        if(k<1) throw std::domain_error("Bad tolerance");
        if(eps<std::fabs(c[k])) {m=k; break;}
    }
}

double ChebyshevApproximation:: operator()(double x) const
{
    if(x<xmin || x>xmax) throw std::domain_error("Bad argument");
    double t=(2*x-xmin-xmax)/(xmax-xmin),p=1,q=t,s=c[0]/2+c[1]*t;
    for(int k=2;k<=m;k++)
    {
        double r=2*t*q-p;
        s+=c[k]*r;
        p=q;
        q=r;
    }
    return s;

}

double ChebyshevApproximation::derivative(double x) const
{
    if(x<xmin || x>xmax) throw std::domain_error("Bad argument");
    double t=(2*x-xmin-xmax)/(xmax-xmin),p=1,q=4*t,s=c[1]+4*c[2]*t;
    for(int k=3;k<=m;k++)
    {
        double r=k*(2*t*q/(k-1)-p/(k-2));
        s+=c[k]*r;
        p=q;
        q=r;
    }
    return 2*s/(xmax-xmin);
}

ChebyshevApproximation ChebyshevApproximation::derivative() const{
    double mi=4/(xmax-xmin);
    std::vector<double>cprim(c.size());
    cprim[m-1]=mi*m*c[m];
    cprim[m-2]=mi*(m-1)*c[m-1];
    for(int k=m-3;k>=0;k--)
        cprim[k]=cprim[k+2]+mi*(k+1)*c[k+1];
    return ChebyshevApproximation(cprim,xmin,xmax);
}

ChebyshevApproximation ChebyshevApproximation::antiderivative() const{
    std::vector<double>cprim(m+2);
    for(int k=1;k<m;k++)
        cprim[k]=(((xmax-xmin)/4)/k)*(c[k-1]-c[k+1]);
    cprim[m]=((xmax-xmin)/4)/m*c[m-1];
    cprim[m+1]=((xmax-xmin)/4)/(m+1)*c[m];
    return ChebyshevApproximation(cprim,xmin,xmax);
}

double ChebyshevApproximation::integrate(double a, double b) const{
    if(a<xmin || b<xmin || a>xmax || b>xmax) throw std::domain_error("Bad interval");
    ChebyshevApproximation f=this->antiderivative();
    return f(b)-f(a);
}

double ChebyshevApproximation::integrate() const{
//3.str p10
    double s=0;
    for(int k=1;k<=(m+1)/2;k++)
        s+=2*c[2*k]/(1-4*k*k);
    s*=(xmax-xmin)/2;
    s+=((xmax-xmin)/2)*c[0];
    return s;

}



int main()
{
    const double PI12=4*std::atan(1);
    auto funsin2=[](double x) {return std::sin(x);};
    ChebyshevApproximation sinch2(funsin2,0,PI12,10);
    std::cout<<sinch2.integrate(0,PI12/2)<<" "<<sinch2.integrate();

    std::cout<<std::endl;

    try
    {
        ChebyshevApproximation cheby([](double x){return x;},1,0,10);
    }
    catch(std::domain_error e)
    {
        std::cout<<e.what()<<std::endl;
    }

    try
    {
        ChebyshevApproximation cheby([](double x){return x;},5,0,0);
    }
    catch(std::domain_error e)
    {
        std::cout<<e.what()<<std::endl;
    }

    try
    {
        ChebyshevApproximation cheby([](double x){return x;},1,5,4);
        cheby.set_m(1);
    }
    catch(std::domain_error e)
    {
        std::cout<<e.what()<<std::endl;
    }

    try
    {
        ChebyshevApproximation cheby([](double x){return x;},1,5,4);
        cheby.trunc(-1);
    }
    catch(std::domain_error e)
    {
        std::cout<<e.what()<<std::endl;
    }

    try
    {
        ChebyshevApproximation cheby([](double x){return x;},1,5,4);
        cheby.derivative(-1);
    }
    catch(std::domain_error e)
    {
        std::cout<<e.what()<<std::endl;
    }

    try
    {
        ChebyshevApproximation cheby([](double x){return x;},1,5,4);
        cheby.integrate(-1,-2);
    }
    catch(std::domain_error e)
    {
        std::cout<<e.what()<<std::endl;
    }



    return 0;
}