#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <iomanip>

class Vector{
        public:
        std::vector<double>elements;
 
        explicit Vector(int n);
        Vector(std::initializer_list<double> l);
        int NElems() const;
        double &operator[](int i);
        double operator[](int i) const;
        double &operator()(int i);
        double operator()(int i) const;
        double Norm() const;
        friend double VectorNorm(const Vector &v);
        double GetEpsilon() const;
        void Print(char separator = '\n', double eps = -1) const;
        friend void PrintVector(const Vector &v, char separator ,double eps );
        friend Vector operator +(const Vector &v1, const Vector &v2);
        Vector &operator +=(const Vector &v);
        friend Vector operator -(const Vector &v1, const Vector &v2);
        Vector &operator -=(const Vector &v);
        friend Vector operator *(double s, const Vector &v);
        friend Vector operator *(const Vector &v, double s);
        Vector &operator *=(double s);
        friend double operator *(const Vector &v1, const Vector &v2);
        friend Vector operator /(const Vector &v, double s);
        Vector &operator /=(double s);

};

Vector::Vector(int n)
{
    if (n<=0) throw std::range_error("Bad dimension");
    Vector::elements.resize(n);
    for(int i:elements) i=0;
}

Vector::Vector(std::initializer_list<double> l){
    if (l.size()<=0) throw std::range_error("Bad dimension");
    for(double i:l) elements.push_back(i);
    
}

int Vector:: NElems() const{ return elements.size();}

double &Vector::operator[](int i){
    return elements[i];
}

double Vector:: operator[](int i) const{
    return elements[i];
}

double &Vector:: operator()(int i) {
    if(i-1<0 || i-1>=NElems()) throw std::range_error("Invalid index");
    return elements[i-1];
}

double Vector:: operator()(int i) const{
    if(i-1<0 || i-1>=NElems()) throw std::range_error("Invalid index");
    return elements[i-1];
}

double Vector:: Norm() const{
    double sum=0;
    for(double i:elements) sum+=i*i;
    return std::sqrt(sum);
}

double VectorNorm(const Vector &v){
    double sum=0;
    for(double i:v.elements) sum+=i*i;
    return std::sqrt(sum);
}

double Vector:: GetEpsilon() const{
    double eps=std::numeric_limits<double>::epsilon();
    return 10*Norm()*eps;
}

void Vector::Print(char separator, double eps) const
{
    double prag=0;
    if(eps<0) prag=GetEpsilon();

    for(int i=0;i<NElems();i++) {
        if (i==NElems()-1 && elements[i]!=separator) {
            if(std::abs(elements[i])<prag) {
                std::cout<<'0'<<separator; break;
            }
            std::cout<<elements[i]; break;
        }
        if(std::abs(elements[i])<eps) {
            std::cout<<'0'<<separator; continue;
        }
        std::cout<<elements[i]<<separator;
    }
}

void PrintVector(const Vector &v, char separator = '\n', double eps = -1){
    double prag=0;
    if(eps<0) prag=v.GetEpsilon();

    for(int i=0;i<v.NElems();i++) {
        if (i==v.NElems()-1 && v.elements[i]!=separator) {
            if(std::abs(v.elements[i])<prag) {
                std::cout<<'0'<<separator; break;
            }
            std::cout<<v.elements[i]; break;
        }
        if(std::abs(v.elements[i])<eps) {
            std::cout<<'0'<<separator; continue;
        }
        std::cout<<v.elements[i]<<separator;
    }
}

Vector operator +(const Vector &v1, const Vector &v2){
    if(v1.NElems()!=v2.NElems()) throw std:: domain_error ("Incompatible formats");

    Vector v3(v1.NElems());
    for(int i=0;i<v1.NElems();i++)
        v3.elements[i]=v1.elements[i]+v2.elements[i];

    return v3;
}

Vector &Vector:: operator +=(const Vector &v)
{
    if(v.NElems()!=NElems()) throw std:: domain_error ("Incompatible formats");

    Vector sum=*this;
    for(int i=0;i<NElems();i++)
        sum.elements[i]=sum.elements[i]+v.elements[i];
    *this=sum;

    return *this;
}

Vector operator -(const Vector &v1, const Vector &v2){
    if(v1.NElems()!=v2.NElems()) throw std:: domain_error ("Incompatible formats");

    Vector v3(v1.NElems());
    for(int i=0;i<v1.NElems();i++)
        v3.elements[i]=v1.elements[i]-v2.elements[i];
        
    return v3;
}

Vector &Vector:: operator -=(const Vector &v)
{
    if(v.elements.size()!=elements.size()) throw std:: domain_error ("Incompatible formats");

    Vector sum=*this;
    for(int i=0;i<v.NElems();i++)
        sum.elements[i]=sum.elements[i]-v.elements[i];

    *this=sum;
    return *this;
}

Vector operator *(double s, const Vector &v){
    Vector res=v;
    for(int i=0;i<v.NElems();i++)
        res.elements[i]=s*res.elements[i];
    return res;
}

Vector operator *(const Vector &v, double s){
    Vector res=v;
    for(int i=0;i<v.NElems();i++)
        res.elements[i]*=s;
    return res;
}

Vector &Vector:: operator *=(double s)
{
    Vector sum=*this;
    for(int i=0;i<NElems();i++)
        sum.elements[i]=sum.elements[i]*s;
    *this=sum;
    return *this;
}

double operator *(const Vector &v1, const Vector &v2){
    if(v1.NElems()!=v2.NElems()) throw std:: domain_error ("Incompatible formats");

    double res=0;
    for(int i=0;i<v1.NElems();i++)
        res+=v1.elements[i]*v2.elements[i];
    return res;

}

Vector operator /(const Vector &v, double s){
    if(s==0) throw std::domain_error ("Division by zero");

    Vector res(v.NElems());
    for(int i=0;i<v.NElems();i++)
        res.elements[i]=v.elements[i]/s;

    return res;
}

Vector &Vector:: operator /=(double s)
{
    if(s==0) throw std::domain_error ("Division by zero");

    Vector sum=*this;
    for(int i=0;i<NElems();i++)
        sum.elements[i]=sum.elements[i]/s;
    
    *this=sum;
    return *this;
}



class Matrix{
    private:
        std::vector<std::vector<double>>elements;
    public:
        Matrix(int m, int n);
        Matrix(const Vector &v);
        Matrix(std::initializer_list<std::vector<double>> l);
        int NRows() const;
        int NCols() const;
        double *operator[](int i);
        const double *operator[](int i) const;
        double &operator()(int i, int j);
        double operator()(int i, int j) const;
        double Norm() const;
        friend double MatrixNorm(const Matrix &m);
        double GetEpsilon() const;
        void Print(int width = 10, double eps = -1) const;
        friend void PrintMatrix(const Matrix &m, int width , double eps );
        friend Matrix operator +(const Matrix &m1, const Matrix &m2);
        Matrix &operator +=(const Matrix &m);
        friend Matrix operator -(const Matrix &m1, const Matrix &m2);
        Matrix &operator -=(const Matrix &m);
        friend Matrix operator *(double s, const Matrix &m);
        friend Matrix operator *(const Matrix &m, double s);
        Matrix &operator *=(double s);
        friend Matrix operator *(const Matrix &m1, const Matrix &m2);
        Matrix &operator *=(const Matrix &m);
        friend Vector operator *(const Matrix &m, const Vector &v);
        friend Matrix Transpose(const Matrix &m);
        void Transpose();

};

Matrix::Matrix(int m,int n){
    if (m<=0||n<=0) throw std::range_error("Bad dimension");
    Matrix::elements.resize(m,std::vector<double>(n));
    for(int i=0;i<m;i++)
        for(int j=0;j<n;j++)
            elements[i][j]=0;
}

Matrix::Matrix(const Vector &v){
    for(int i=0;i<v.NElems();i++)
        elements[i][0]=v.elements[i];
}

Matrix::Matrix(std::initializer_list<std::vector<double>> l){
    if (l.size()<=0 || l.begin()->size()<=0) throw std::range_error("Bad dimension");
    
    for(auto i=l.begin();i<l.end();i++) if(i->size()!=l.begin()->size()) throw std::logic_error("Bad matrix");

    for (const auto &i : l) elements.push_back(i);
}

int Matrix::NRows() const{
    return elements.size();
}

int Matrix::NCols() const{
    if(elements.size() != 0) return elements[0].size();
    return 0;
}

double *Matrix::operator[](int i){
    std::vector<double>&red=elements[i];
    double *r=&red[0];
    return r;
}

const double *Matrix:: operator[](int i) const{
    const std::vector<double>&red=elements[i];
    const double *r=&red[0];
    return r;
}

double &Matrix::operator()(int i, int j){
    if (i-1<0 || i-1>=NRows() || j-1<0 || j-1>=NCols()) throw std::range_error("Invalid index");
    return elements[i-1][j-1];
}

double Matrix::operator()(int i, int j) const{
    if (i-1<0 || i-1>=NRows() || j-1<0 || j-1>=NCols()) throw std::range_error("Invalid index");
    return elements[i-1][j-1];
}

double Matrix::Norm() const{
    double sum=0;
    for (int i=0;i<NRows();i++)
        for (int j=0;j<NCols();j++)
            sum+=elements[i][j]*elements[i][j];
    
    double res=sqrt(sum);
    return res; 
}

double MatrixNorm(const Matrix &m){
    double sum=0;
    for (int i=0;i<m.NRows();i++)
        for (int j=0;j<m.NCols();j++)
            sum+=m.elements[i][j]*m.elements[i][j];
    
    double res=sqrt(sum);
    return res; 
}

double Matrix::GetEpsilon() const
{
    double eps=std::numeric_limits<double>::epsilon();
    return 10*Norm()*eps;
}

void Matrix::Print(int width, double eps) const{
    double prag=0;
    if(eps<0) prag=GetEpsilon();
    std::cout<<std::setw(width);
    for(int i=0;i<NRows();i++) {
        for(int j=0;j<NCols();j++)
        {
            if(std::abs(elements[i][j])<eps) {
                std::cout<<'0'<<std::setw(width); continue;
            }
            std::cout<<elements[i][j]<<std::setw(width);
        }
        std::cout<<std::endl;
    }
}

void PrintMatrix(const Matrix &m, int width = 10, double eps = -1){
    double prag;
    if(eps<0) prag=m.GetEpsilon();
    std::cout<<std::setw(width);
    for(int i=0;i<m.NRows();i++) {
        for(int j=0;j<m.NCols();j++)
        {
            if(std::abs(m.elements[i][j])<eps) {
                std::cout<<'0'<<std::setw(width); continue;
            }
            std::cout<<m.elements[i][j]<<std::setw(width);
        }
        std::cout<<std::endl;
    }
}

Matrix operator +(const Matrix &m1, const Matrix &m2){
    if((m1.NRows()!=m2.NRows()) || (m1.NCols()!=m2.NCols()))
        throw std:: domain_error ("Incompatible formats");
        
    Matrix res(m1.NRows(),m1.NCols());
        for (int i=0;i<m1.NRows();i++)
            for (int j=0;j<m2.NCols();j++)
                res.elements[i][j]=m1.elements[i][j]+m2.elements[i][j];

    return res;
}

Matrix &Matrix::operator +=(const Matrix &m){
    if((m.NRows()!=NRows()) || (m.NCols()!=NCols()))
        throw std:: domain_error ("Incompatible formats");
    Matrix sum=*this;
    for(int i=0;i<sum.NRows();i++)
        for (int j=0;j<sum.NCols();j++)
            sum.elements[i][j]=sum.elements[i][j]+m.elements[i][j];
    
    *this=sum;
    return *this;
}

Matrix operator -(const Matrix &m1, const Matrix &m2){
    if((m1.NRows()!=m2.NRows()) || (m1.NCols()!=m2.NCols()))
        throw std:: domain_error ("Incompatible formats");
        
    Matrix res(m1.NRows(),m1.NCols());
        for (int i=0;i<m1.NRows();i++)
            for (int j=0;j<m2.NCols();j++)
                res.elements[i][j]=m1.elements[i][j]-m2.elements[i][j];

    return res;
}

Matrix &Matrix::operator -=(const Matrix &m){
    if((m.NRows()!=NRows()) || (m.NCols()!=NCols()))
        throw std:: domain_error ("Incompatible formats");

    Matrix sum=*this;
    for(int i=0;i<NRows();i++)
        for (int j=0;j<NCols();j++)
            sum.elements[i][j]=sum.elements[i][j]-m.elements[i][j];
    
    *this=sum;
    return *this;
}

Matrix operator *(double s, const Matrix &m)
{
    Matrix res(m.NRows(),m.NCols());
        for (int i=0;i<m.NRows();i++)
            for (int j=0;j<m.NCols();j++)
                res.elements[i][j]=s*m.elements[i][j];

    return res;
}

Matrix operator *(const Matrix &m, double s){
    Matrix res(m.NRows(),m.NCols());
        for (int i=0;i<m.NRows();i++)
            for (int j=0;j<m.NCols();j++)
                res.elements[i][j]=m.elements[i][j]*s;

    return res;
}

Matrix &Matrix::operator *=(double s){
    Matrix sum=*this;
    for(int i=0;i<sum.NRows();i++)
        for (int j=0;j<sum.NCols();j++)
            sum.elements[i][j]=sum.elements[i][j]*s;
    
    *this=sum;
    return *this;
}

Matrix operator *(const Matrix &m1, const Matrix &m2){
    if((m1.NCols()!=m2.NRows()))
        throw std:: domain_error ("Incompatible formats");
        
    Matrix res(m1.NRows(),m2.NCols());
        for (int i=0;i<m1.NRows();i++)
            for (int j=0;j<m2.NCols();j++)
                    for (int k=0;k<m2.NCols();k++)
                        res.elements[i][j]+=m1.elements[i][k]*m2.elements[k][j];

    return res;
}

Matrix &Matrix:: operator *=(const Matrix &m){
    if((m.NRows()!=NCols()))
        throw std:: domain_error ("Incompatible formats");

    Matrix res(NRows(),m.NCols());

    for (int i=0;i<NRows();i++)
        for (int j=0;j<m.NCols();j++)
                for (int k=0;k<m.NCols();k++)
                    res.elements[i][j]+=elements[i][k]*m.elements[k][j];

    *this=res;
    return *this;
}

Vector operator *(const Matrix &m, const Vector &v){
    if((m.NCols()!=v.NElems()))
        throw std:: domain_error ("Incompatible formats");

    Vector res(m.NRows());

    for(int i=0;i<m.NRows();i++)
        for(int j=0;j<v.NElems();j++)
            res.elements[i]+=m.elements[i][j]*v.elements[j];
            
    return res;
}

Matrix Transpose(const Matrix &m){
    Matrix transp(m.NRows(),m.NCols());

    if (m.NRows()!=m.NCols()){
        Matrix transp(m.NCols(),m.NRows());
        for(int i=0;i<m.NRows();i++)
            for(int j=0;j<m.NCols();j++)
                transp[j][i]=m.elements[i][j];
        return transp;

    } 

    for(int i=0;i<m.NRows();i++)
        for(int j=0;j<m.NCols();j++)
            {
                double temp=m.elements[j][i]; 
                transp.elements[j][i]=m.elements[i][j]; 
                transp.elements[i][j]=temp;
            }

    return transp;
}

void Matrix::Transpose()
{
    if (NRows()!=NCols()){
        Matrix transp(NCols(),NRows());
        for(int i=0;i<NRows();i++)
            for(int j=0;j<NCols();j++)
                transp[j][i]=elements[i][j];
        *this=transp;
    }
    else {
        for(int i=0;i<NRows();i++)
            for(int j=i+1;j<NCols();j++)
                std::swap(elements[i][j],elements[j][i]);
    }
}

int main()
{
    Vector v1(2);
    v1={2,3};
    std::cout<<"Vector v1(2):"<<std::endl;
    std::cout<<std::endl;
    std::cout<<"Sekvencijski konstruktor Print() ";
    v1.Print(); std::cout<<std::endl;
    std::cout<<"Sekvencijski konstruktor PrintVector(v1): "<<std::endl;
    PrintVector(v1);
    std::cout<<"NElems(): "<<v1.NElems()<<std::endl;
    std::cout<<"operator[0]: "<<v1[0]<<std::endl;
    std::cout<<"&operator[1]: "<<v1[1]<<std::endl;
    std::cout<<"operator(1):" <<v1(1)<<std::endl;
    std::cout<<"&operator(2): "<<v1(2)<<std::endl;
    std::cout<<"Norm(): "<<v1.Norm()<<std::endl;
    std::cout<<"Vector Norm(v1): "<<VectorNorm(v1)<<std::endl;
    std::cout<<"GetEpsilon(): "<<v1.GetEpsilon()<<std::endl<<std::endl<<std::endl;

    const Vector v2=v1;
    std::cout<<"Vector v2(2):"<<std::endl;
    std::cout<<std::endl;
    std::cout<<"Sekvencijski konstruktor v2.Print(' ') ";
    v2.Print(' '); std::cout<<std::endl;
    std::cout<<"Sekvencijski konstruktor PrintVector(v2,' ',10): "<<std::endl;
    PrintVector(v2,' ',10);
    std::cout<<"NElems(): "<<v2.NElems()<<std::endl;
    std::cout<<"operator[0]: "<<v2[0]<<std::endl;
    std::cout<<"&operator[1]: "<<v2[1]<<std::endl;
    std::cout<<"operator(1):" <<v2(1)<<std::endl;
    std::cout<<"&operator(2): "<<v2(2)<<std::endl;
    std::cout<<"Norm(): "<<v2.Norm()<<std::endl;
    std::cout<<"Vector Norm(v2): "<<VectorNorm(v2)<<std::endl;
    std::cout<<"GetEpsilon(): "<<v2.GetEpsilon()<<std::endl;
    std::cout<<"operator +(v1,v2): "<<std::endl;
    Vector sum=v1+v2; sum.Print();
    std::cout<<"operator -(v1,v2): "<<std::endl;
    sum=v1-v2; sum.Print();
    std::cout<<"operator /(v1,2): "<<std::endl;
    sum=v1/2; sum.Print();
    std::cout<<"operator *(v1,2): "<<std::endl;
    sum=v1*2; sum.Print();
    std::cout<<"operator *(3,v1): "<<std::endl;
    sum=3*v1; sum.Print();
    std::cout<<"operator *(v1,v2): "<<v1*v2<<std::endl;
    std::cout<<"operator +=(v2): "<<std::endl;
    v1+=v2; v1.Print();
    std::cout<<"operator -=(v2): "<<std::endl;
    v1-=v2;v1.Print();
    std::cout<<"operator *=(10): "<<std::endl;
    v1*=10;v1.Print();
    std::cout<<"operator /=(10): "<<std::endl<<std::endl<<std::endl;
    v1/=10;v1.Print();

    std::cout<<"Matrica m1(3,3):"<<std::endl;
    Matrix m1(3,3);
    m1={{1,1,1},{2,2,2},{3,3,3}};
    std::cout<<std::endl;
    std::cout<<"Sekvencijski konstruktor m1.Print(): ";
    m1.Print(); std::cout<<std::endl;
    std::cout<<"Sekvencijski konstruktor PrintMatrix(m1): ";
    PrintMatrix(m1); std::cout<<std::endl;
    std::cout<<"Sekvencijski konstruktor m1.Print(4): ";
    m1.Print(4); std::cout<<std::endl;
    std::cout<<"Sekvencijski konstruktor PrintMatrix(m1,4): ";
    PrintMatrix(m1,4); std::cout<<std::endl;
    const Matrix m2=m1;
    std::cout<<"Sekvencijski konstruktor PrintMatrix(m2,4,50): ";
    PrintMatrix(m2,4,50); std::cout<<std::endl;
    std::cout<<"NRows(): "<<m2.NRows()<<std::endl;
    std::cout<<"NCols(): "<<m2.NCols()<<std::endl;
    std::cout<<"Norm() "<<m2.Norm()<<std::endl;
    std::cout<<"operator(1,2):" <<m2(1,2)<<std::endl;
    std::cout<<"&operator(1,2): "<<m2(1,2)<<std::endl;
    std::cout<<"operator[1][2]: " <<m2[0][1]<<std::endl;
    std::cout<<"&operator[1][2]: "<<m2[0][1]<<std::endl;
    std::cout<<"Matrix Norm: "<<MatrixNorm(m2)<<std::endl;
    std::cout<<"GetEpsilon: "<<m2.GetEpsilon()<<std::endl;
    std::cout<<"Matrix: "<<m2.GetEpsilon()<<std::endl;
    std::cout<<"operator +(m1,m2): "<<std::endl;
    Matrix res=m1+m2; res.Print();
    std::cout<<"operator -(m1,m2): "<<std::endl;
    res=m1-m2; res.Print();
    std::cout<<"operator *(m1,2): "<<std::endl;
    res=m1*2; res.Print();
    std::cout<<"operator *(3,m1): "<<std::endl;
    res=3*m1; res.Print();
    std::cout<<"operator *(m1,m2): "<<std::endl;
    res=m1*m2; res.Print();
    Vector v5(3); v5={2,2,2};
    std::cout<<"operator *(m1,v5): "<<std::endl;
    Vector res2(3);
    res2=m1*v5; res2.Print();
    std::cout<<"operator +=(m2): "<<std::endl;
    m1+=m2; m1.Print();
    std::cout<<"operator -=(m2): "<<std::endl;
    m1-=m2;m1.Print();
    std::cout<<"operator *=(10): "<<std::endl;
    m1*=10;m1.Print();
    std::cout<<"operator *=(m2): "<<std::endl;
    m1*=m2;m1.Print();
    std::cout<<"m1.Transpose(): "<<std::endl;
    m1.Transpose();
    PrintMatrix(m1);
    std::cout<<"Transpose(m2): "<<std::endl<<std::endl;
    Transpose(m2);
    PrintMatrix(m2);
    std::cout<<"Matrix m4: "<<std::endl;
    Matrix m4(5,3);
    m4={{1,2,3},{4,5,6},{7,8,9},{10,11,12},{13,14,15}};
    m4.Print();
    std::cout<<"Transpose(m4): "<<std::endl;
    m4.Transpose();
    PrintMatrix(m4);
   

//Izuzeci:

    // Vector test(0); 
    // Vector test(-1);
    // Vector test={{}};
    // Vector test={};
    // Vector test1={1};
    // test1(2);
    // test1(0);
    // test1/0;
    // Vector test2={2,3};
    // test1*test2;
    // test1+test2;
    // test1-test2;

    // Matrix proba(0,0);
    // Matrix proba(-1,-1);
    // Matrix proba(-1,2);
    // Matrix proba(2,-1);
    // Matrix proba={{2,3},{1,1,1},{3,3,3}}; // Grabava matrica
    // Matrix proba={};
    // Matrix proba={{}};
    // Matrix proba ={{1},{1}};
    // proba(0,0);
    // proba(2,0);
    // Matrix proba2={{1,2},{3,4}};
    // proba*proba2;
    // proba-proba2;
    // proba+proba2;

    return 0;
}


