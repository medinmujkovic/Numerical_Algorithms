#include <iostream>
#include <stdexcept>
#include <vector>
#include <algorithm>
#include <cmath>

class AbstractInterpolator{
    private:
        mutable std::vector<std::pair<double,double>>data;
        mutable int k=1;
    protected:
        int Locate(double x) const;
    public:
        AbstractInterpolator(const std::vector<std::pair<double, double>> &data);
        virtual double operator()(double x) const = 0;
        const std::vector<std::pair<double, double>> &getData() const {return data;}
    
};

    AbstractInterpolator::AbstractInterpolator(const std::vector<std::pair<double, double>> &data)
    {
        this->data=data;
        std::sort(this->data.begin(),this->data.end(),[](const std::pair<double,double> &a, const std::pair<double,double> &b)
        {
            //if(a.first==b.first) throw std::domain_error("Invalid data set");
            return a.first<b.first;
        });

        for (int i=0;i<data.size()-1;i++)
            if(fabs(data[i].first-data[i+1].first)<std::numeric_limits<double>::epsilon()) throw std::domain_error("Invalid data set");       
    }

    int AbstractInterpolator::Locate(double x) const
    {      
        if (x<=data[0].first) return 0;
        if (x>data[data.size()-1].first) return data.size();
        if(data[k-1].first<=x && data[k].first>x) return k;
        std::pair<double,double>temp(x,0);
        auto low=std::lower_bound(data.begin(),data.end(),temp,[](const std::pair<double,double> &a, const std::pair<double,double> &b){
             return a.first<b.first;
        });
        k=low-data.begin();
        return k;
    }

class LinearInterpolator : public AbstractInterpolator{
    public:
    LinearInterpolator(const std::vector<std::pair<double, double>> &data):AbstractInterpolator(data){}
    double operator()(double x) const override;
};
    double Pomocna(double x,std::pair<double,double>t1,std::pair<double,double>t2)
    {
        return t1.second*((t2.first-x)/(t2.first-t1.first))+t2.second*((x-t1.first)/(t2.first-t1.first));
    }
    double LinearInterpolator:: operator()(double x) const
    {
        int index=Locate(x);
        auto data=getData();
        if(index<=1) return Pomocna(x,data[0],data[1]);
        else if(index>=data.size()) return Pomocna(x,data[data.size()-2],data[data.size()-1]);
        return Pomocna(x,data[index-1],data[index]);
    } 

class PolynomialInterpolator : public AbstractInterpolator{
    private:
        std::vector<double>q;
    public:
    PolynomialInterpolator(const std::vector<std::pair<double, double>> &data);
    double operator()(double x) const override;
    void AddPoint(const std::pair<double, double> &p);
    std::vector<double> GetCoefficients() const;
};

    PolynomialInterpolator::PolynomialInterpolator(const std::vector<std::pair<double, double>> &data):AbstractInterpolator(data)
    {
        for (int i=0;i<data.size();i++)
            q.push_back(data[i].second);

        for (int j=0;j<data.size();j++)
            for (int i=data.size()-1;i>=j+1;i--){
                q[i]=(q[i]-q[i-1])/(data[i].first-data[i-j-1].first);
            } 
    }

    double PolynomialInterpolator::operator()(double x) const
    {
        auto data=getData();
        double result=q[q.size()-1];

        for (int i=q.size()-2;i>=0;i--) 
            result=result*(x-data[i].first)+q[i];

        return result;
    }
    
    void PolynomialInterpolator::AddPoint(const std::pair<double, double> &p) {
        auto data=getData();
        int n=data.size();
        if (std::find_if(data.begin(), data.end(), [p](const std::pair<double, double> &a) {
                return a.first==p.first;
            })!=data.end()) {
                throw std::domain_error("Invalid point");
            }
        data.push_back(p);
        q.resize(n+1);
        double temp=0;
        for (int i=1;i<n;i++) {
            temp=q[i-1];
            q[n-1]=data[n-1].second;
            data[n-1].second=(data[i-1].second-temp)/(data[n-1].first-data[n-1-i].first);
        }
            q[n-1]=data[n-1].second;
    }

    std::vector<double> PolynomialInterpolator::GetCoefficients() const{
        auto data=getData();
        int n=data.size();
        std::vector<double>w(n+1);
        std::vector<double>p(n+1);
        p[0]=q[0];
        w[0]=1;
        
        for (int i=1;i<=n;i++) {
            w[i]=w[i-1];
            for (int j=i-1;j>=1;j--)
                w[j]=w[j-1]-data[i-1].first*w[j];
            w[0]=-data[i-1].first*w[0];
            for (int j=0;j<n;j++)
                p[j]+=q[i]*w[j];
        }

        return p;
    }

class PiecewisePolynomialInterpolator : public AbstractInterpolator{
    private:
        int order;
    public:
        PiecewisePolynomialInterpolator(const std::vector<std::pair<double, double>> &data,int order);
        double operator()(double x) const override;
};

   PiecewisePolynomialInterpolator::PiecewisePolynomialInterpolator(const std::vector<std::pair<double, double>> &data,int order):AbstractInterpolator(data),
        order(order)
    {
        if (order<1||order>=data.size()) throw std::domain_error("Invalid order");

    }
    
    double PiecewisePolynomialInterpolator::operator()(double x) const{
        int index=Locate(x);
        double x0,y0;
        auto data=getData();
        double result=0;

        if (order%2==0) {
            x0=index-order/2-1;
            y0=index+order/2;
        } 
        else {
            x0=index-(order-1)/2-1;
            y0=index+(order+1)/2;
        }

        if (y0>=data.size())
        {
            x0=data.size()-order-1;
            y0=data.size();
        }
        if(x0<=0)
        {
            x0=0;
            y0=order+1;
        }

        result=0;
        for (int i=x0;i<y0;i++)
        {
            double temp=data[i].second;

            for (int j=x0;j<y0;j++)
                if (j!=i)
                    temp=temp*(x-data[j].first)/(data[i].first-data[j].first);
            result+=temp;
        }
    
        return result;   
    }

class SplineInterpolator : public AbstractInterpolator{
    private:
        std::vector<double>alfa;
        std::vector<double>r;
        std::vector<double>q;
    public:
        SplineInterpolator(const std::vector<std::pair<double, double>> &data);
        double operator()(double x) const override;
};

    SplineInterpolator::SplineInterpolator(const std::vector<std::pair<double, double>> &data) : AbstractInterpolator(data) {
        int n=data.size();
        r.resize(n); r[0]=0; r[n-1]=0;
        alfa.resize(n);
        q.resize(n);

        for(int i=1;i<n-1;i++)
        {
            alfa[i]=2*(data[i+1].first-data[i-1].first);
            r[i]=3*((data[i+1].second-data[i].second)/(data[i+1].first-data[i].first)-(data[i].second-data[i-1].second)/(data[i].first-data[i-1].first));
        }
        for(int i=1;i<n-2;i++)
        {
            double mi=(data[i+1].first-data[i].first)/alfa[i];
            alfa[i+1]-=mi*(data[i+1].first-data[i].first);
            r[i+1]-=mi*r[i];
        }
        r[n-2]/=alfa[n-2];
        for(int i=n-3;i>0;i--)
            r[i]=(r[i]-r[i+1]*(data[i+1].first-data[i].first))/alfa[i];
        
        for(int i=0;i<n-1;i++)
        {
            double delta=data[i+1].first-data[i].first;
            alfa[i]=(r[i+1]-r[i])/(3*delta);
            q[i]=(data[i+1].second-data[i].second)/delta-delta*(r[i+1]+2*r[i])/3;
        }

        
    }

double SplineInterpolator::operator()(double x) const {
    int index=Locate(x);
    auto data=getData();
    int n=data.size();

    if (index<=0) {
        index=1;
    } else if (index>=n) {
        index=n-1;
    }

    double t=x-data[index-1].first;
    return data[index-1].second+t*(q[index-1]+t*(r[index-1]+t*alfa[index-1]));
}




class BarycentricInterpolator : public AbstractInterpolator{
    private:
        int order;
        std::vector<double>w;
    public:
        BarycentricInterpolator(const std::vector<std::pair<double, double>> &data, int order);
        double operator()(double x) const override;
        std::vector<double> GetWeights() const{return w;}
};

    BarycentricInterpolator::BarycentricInterpolator(const std::vector<std::pair<double, double>> &data, int order):AbstractInterpolator(data){
        if (order<0 || order>=data.size()) throw std::domain_error("Invalid order");
        int n=data.size();
        w.resize(n,0);
        for (int i=0;i<n;i++)
        {
            w[i]=0;
            for (int k=std::max(0,i-order);k<=std::min(i,n-order-1);k++)
            {
                double p=1;
                for(int j=k;j<(k+order);j++)
                {
                    if(j!=i) p=p/(data[i].first-data[j].first);
                }
                if (k%2==1) p=-p;

                w[i]=w[i]+p;
            }
            
        }
        
    }

    double BarycentricInterpolator::operator()(double x) const{
        double p=0;
        double q=0;
        auto data=getData();
        int n=data.size();

        for(int i=0;i<n;i++)
        {
            if(x==data[i].first) return data[i].second;
            double u=w[i]/(x-data[i].first);
            p+=u*data[i].second;
            q+=u;
        }
        return p/q;
    }


class TrigonometricInterpolator : public AbstractInterpolator{
    public:
        TrigonometricInterpolator(const std::vector<std::pair<double, double>> &data);
        double operator()(double x) const override;
};

    TrigonometricInterpolator::TrigonometricInterpolator(const std::vector<std::pair<double, double>> &data):AbstractInterpolator(data)
    {
        if(data.front().first!=data.back().first) throw std::domain_error("Function is not periodic");
    }
    double TrigonometricInterpolator::operator()(double x) const{
        auto data=getData();
        int n=data.size();
        double res=0;
        for (int k=1;k<=n-1;k++)
        {
            double pr=data[k].second;
            double c=1;
            for(int j=1;j<=n-1;j++)
            {
                double a=0;
                double b=0;
                if(n%2==0 && k!=j)
                {
                    double a=sin(4*atan(1)/n)*(x-data[j].first);
                    double b=sin(4*atan(1)/n)*(data[k].first-data[j].first);
                    c=a/b; pr*=c;
                }
                else if(n%2!=0 && k!=j)
                {
                    std::vector<double>a(n,0);
                    std::vector<double>b(n,0);
                    c=a[0]/2+a[n]*cos((2*4*atan(1))/x);
                    for(int i=1;i<n-1;i++)
                        c+=a[k]*cos((2*4*atan(1))/x)+b[k]+sin((2*4*atan(1))/x);
                    pr*=c;
                }
                
            }
            res+=data[k].second*pr;
        }
        return res;
    }

int main()
{
    const double PI4=std::atan(1)*4;
    std::vector<std::pair<double,double>> data4;
    for(double i=2*PI4; i>=0; i-=PI4/2)
        data4.push_back({i,std::cos(i)});
    SplineInterpolator si4(data4);
    std::cout<<si4(-0.1)<<" "<<std::cos(-0.1)<<std::endl;
    std::cout<<si4(-0.2)<<" "<<std::cos(-0.2)<<std::endl;
    std::cout<<std::round(si4(PI4/2))<<" "<<std::round(std::cos(PI4/2))<<std::endl;
    std::cout<<si4(PI4/2+0.1)<<" "<<std::cos(PI4/2+0.1)<<std::endl;
    std::cout<<si4(PI4*3+0.1)<<" "<<std::cos(PI4*3+0.1)<<std::endl;
    std::cout<<si4(PI4*3+0.2)<<" "<<std::cos(PI4*3+0.2);
    std::vector<std::pair<double, double>> data = {
        {1.2, 3.4},
        {5.6, 7.8},
        {9.0, 10.1},
        {11.2, 13.4},
        {15.6, 11.8}
    };
    LinearInterpolator li({{1,4},{3,6},{7,8}});
    PolynomialInterpolator pi({{1,1},{2,2},{3,3}});
    std::cout<<li(2)<<" "<<li(5.5)<<" "<<li(6.25);

    return 0;
}

