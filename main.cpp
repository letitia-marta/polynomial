#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <iomanip>
#define show(f) std::cout<<"Fie polinomul "<<#f<<". "<<#f<<"(X) = ";
using namespace std;
int cmmdc (int a, int b)
{
    while (b)
    {
        int r=a%b;
        a=b;
        b=r;
    }
    return a;
}

class NrRational
{
    int numarator,numitor;
    public:
        NrRational();
        NrRational(int,int);
        int getNumarator()
        {
            return numarator;
        }
        int getNumitor()
        {
            return numitor;
        }
        float fractie=1.0*numarator/numitor;
        bool operator > (float k);
        bool operator == (float k);
        bool operator < (float k);
        NrRational& operator = (int x);
        friend NrRational operator + (const NrRational&, const NrRational&);
        friend NrRational operator * (const NrRational&, const NrRational&);
        friend istream& operator >> (istream&, NrRational&);
        friend ostream& operator << (ostream&, NrRational&);
};
NrRational::NrRational()
{
    this->numarator=0;
    this->numitor=1;
}
NrRational::NrRational (int a, int b)
{
    this->numarator=a;
    this->numitor=b;
}
bool NrRational::operator > (float k)
{
    return (this->fractie>k);
}
bool NrRational::operator == (float k)
{
    return (this->fractie==k);
}
bool NrRational::operator < (float k)
{
    return (this->fractie<k);
}
NrRational& NrRational::operator = (int x)
{
    this->numarator=x;
    this->numitor=1;
    return *this;
}
NrRational operator + (const NrRational &a, const NrRational &b)
{
    NrRational s;
    int n1,n2;
    n1=a.numarator*b.numitor+a.numitor*b.numarator;
    n2=a.numitor*b.numitor;
    s.numarator=n1/cmmdc(n1,n2);
    s.numitor=n2/cmmdc(n1,n2);
    return s;
}
NrRational operator * (const NrRational &a, const NrRational &b)
{
    NrRational p;
    int n1,n2;
    n1=a.numarator*b.numarator;
    n2=a.numitor*b.numitor;
    p.numarator=n1/cmmdc(n1,n2);
    p.numitor=n2/cmmdc(n1,n2);
    return p;
}
istream& operator >> (istream& in, NrRational& a)
{
    in>>a.numarator>>a.numitor;
    return in;
}
ostream& operator << (ostream& out, NrRational& a)
{
    out<<a.getNumarator();
    if (a.getNumitor()!=1)
        out<<'/'<<a.getNumitor();
    return out;
}

class NrComplex
{
    int re,im;
    public:
        NrComplex();
        NrComplex(int,int);
        int getRe()
        {
            return re;
        }
        int getIm()
        {
            return im;
        }
        bool operator > (int k);
        bool operator == (int k);
        bool operator < (int k);
        NrComplex& operator = (int x);
        friend NrComplex operator + (const NrComplex&, const NrComplex&);
        friend NrComplex operator * (const NrComplex&, const NrComplex&);
        friend istream& operator >> (istream&, NrComplex&);
        friend ostream& operator << (ostream&, const NrComplex&);
};
NrComplex::NrComplex()
{
    this->re=0;
    this->im=0;
}
NrComplex::NrComplex (int real, int imag)
{
    this->re=real;
    this->im=imag;
}
bool NrComplex::operator > (int k)
{
    return true;
}
bool NrComplex::operator == (int k)
{
    return false;
}
bool NrComplex::operator < (int k)
{
    return (this->re<k);
}
NrComplex& NrComplex::operator = (int x)
{
    this->re=x;
    this->im=0;
    return *this;
}
NrComplex operator + (const NrComplex &a, const NrComplex &b)
{
    NrComplex s;
    s.re=a.re+b.re;
    s.im=a.im+b.im;
    return s;
}
NrComplex operator * (const NrComplex &a, const NrComplex &b)
{
    NrComplex p;
    p.re=a.re*b.re-a.im*b.im;
    p.im=a.im*b.re+a.re*b.im;
    return p;
}
istream& operator >> (istream& in, NrComplex& z)
{
    in>>z.re>>z.im;
    return in;
}
ostream& operator << (ostream& out, const NrComplex& z)
{
    out<<'(';
    if (z.re)
    {
        out<<z.re;
        if (z.im>0)
            out<<'+'<<z.im<<"i)";
        else if (z.im<0)
            out<<z.im<<"i)";
    }
    else
        out<<z.im<<"i)";
    return out;
}

template <typename tip> class Polinom
{
    int grad;
    tip X[100]; //pe pozitia i e coeficientul termenului X^i
    public:
        Polinom();
        void citire();
        void afisare();
        int getGrad() const
        {
            return grad;
        }
        friend Polinom<tip> operator + (Polinom<tip>& p, Polinom<tip>& q)
        {
            Polinom<tip> r;
            r.grad=max(p.grad,q.grad);
            for (int i=0; i<=r.grad; i++)
            {
                int z=0;
                r.X[i]=z;
            }
            for (int i=0; i<=min(p.grad,q.grad); i++)
                r.X[i]=r.X[i]+p.X[i]+q.X[i];
            for (int i=min(p.grad,q.grad)+1; i<=r.grad; i++)
            {
                if (r.grad==p.grad)
                    r.X[i]=r.X[i]+p.X[i];
                else
                    r.X[i]=r.X[i]+q.X[i];
            }
            return r;
        }
        friend Polinom<tip> operator * (Polinom<tip>& p, Polinom<tip>& q)
        {
            Polinom<tip> r;
            r.grad=p.grad+q.grad;
            for (int i=0; i<=r.grad; i++)
                r.X[i]=0;
            for (int i=0; i<=p.grad; i++)
                for (int j=0; j<=q.grad; j++)
                    r.X[i+j]=r.X[i+j]+p.X[i]*q.X[j];
            return r;
        }
        tip operator [] (tip);
        Polinom<tip> derivata();
};
template <typename tip> Polinom<tip>::Polinom()
{
    grad=0;
}
template <typename tip> void Polinom<tip>::citire()
{
    cin>>this->grad;
    for (int i=this->grad; i>=0; i--)
        cin>>this->X[i];
}
template <typename tip> void Polinom<tip>::afisare()
{
    if (this->grad==1)
        cout<<X[this->grad]<<"X";
    else
    {
        if (X[this->grad]>0) //termenul dominant
        {
            if (X[this->grad]==1)
                cout<<"X^"<<this->grad;
            else
                cout<<X[this->grad]<<"X^"<<this->grad;
        }
        else
        {
            if (X[this->grad]==-1)
                cout<<'-'<<"X^"<<this->grad;
            else
                cout<<X[this->grad]<<"X^"<<this->grad;
        }
    }
    for (int i=this->grad-1; i>1; i--)
    {
        if (X[i]>0)
        {
            if (X[i]==1)
                cout<<" +"<<"X^"<<i;
            else
                cout<<" +"<<X[i]<<"X^"<<i;
        }
        else if (X[i]<0)
        {
            if (X[i]==-1)
                cout<<" -"<<"X^"<<i;
            else
                cout<<' '<<X[i]<<"X^"<<i;
        }
    }
    if (X[1]>0)
    {
        if (X[1]==1)
            cout<<" +"<<"X";
        else
            cout<<" +"<<X[1]<<"X";
    }
    else if (X[1]<0)
    {
        if (X[1]==-1)
            cout<<" -"<<"X";
        else
            cout<<' '<<X[1]<<"X";
    }
    if (X[0]>0)
    {
        if (X[0]==1)
            cout<<" +1\n";
        else
            cout<<" +"<<X[0]<<'\n';
    }
    else if (X[0]<0)
    {
        if (X[0]==-1)
            cout<<" -1\n";
        else
            cout<<' '<<X[0]<<'\n';
    }
    else
        cout<<'\n';
}
template <typename tip> tip Polinom<tip>::operator [] (tip x)
{
    tip rez;
    int z=0;
    rez=z;
    for (int i=0; i<=this->grad; i++)
    {
        int u=1;
        tip p;
        p=u;
        for (int j=1; j<=i; j++)
            p=p*x;
        rez=rez+this->X[i]*p;
    }
    return rez;
}
template <typename tip> Polinom<tip> Polinom<tip>::derivata()
{
    Polinom d;
    d.grad=this->grad-1;
    cout<<d.grad<<'\n';
    for (int i=1; i<=this->grad; i++)
    {
        int u=i;
        tip p;
        p=u;
        d.X[i-1]=(this->X[i])*p;
    }
    return d;
}

template <typename T> void alg (Polinom<T>& f, Polinom<T>& g)
{
    cout<<"Dati un polinom:\n";
    f.citire();
    cout<<"Dati alt polinom:\n";
    g.citire();

    Polinom<T> s=f+g;
    show(s);
    s.afisare();
    Polinom<T> p=f*g;
    show(p);
    p.afisare();

    T x;
    cout<<"Dati o valoare: "; cin>>x;
    cout<<"Valoarea polinomului f in "<<x<<" este ";
    T rez;
    rez=f[x];
    cout<<rez;
    cout<<".\nValoarea polinomului g in "<<x<<" este ";
    rez=g[x];
    cout<<rez;
}
void clearscreen()
{
    char ws;
    cout<<".\n\nApasati X pentru a reveni la meniu.\n";
    cin>>ws;
    if (ws=='x' || ws=='X')
        system("CLS");
}

float Newton (Polinom<float> f, float a, float b)
{
    if (f[a]==0)
        return a;
    if (f[b]==0)
        return b;
    float x0,x1;
    Polinom<float> d;
    d=f.derivata();
    x1=f[0]/d[0];
    while (fabs(f[x1]))
    {
        x0=x1;
        x1=x0-(f[x0]/d[x0]);
    }
    return x1;
}
int main()
{
    while (1)
    {
        cout<<"1. Coeficienti intregi\n";
        cout<<"2. Coeficienti reali\n";
        cout<<"3. Coeficienti rationali\n";
        cout<<"4. Coeficienti complecsi\n";
        cout<<"\nIntroduceti varianta dorita, sau apasati X pentru a termina: \n";
        char choice;
        cin>>choice;
        system("CLS");
        if (choice=='1')
        {
            Polinom<int> f,g;
            alg<int>(f,g);
            clearscreen();
        }
        else if (choice=='2')
        {
            Polinom<float> f,g;
            alg<float>(f,g);
            cout<<"\nDati intervalul in care se va cauta radacina lui f(X): ";
            float a,b;
            cin>>a>>b;
            cout<<"Folosim metoda lui Newton.\n";
            cout<<"Calculam succesiv: xn+1 = xn - f(xn)/f'(xn)\n";
            float rad=Newton(f,a,b);
            if (rad>=a && rad<=b)
                cout<<"Radacina este "<<setprecision(1)<<fixed<<rad;
            else
                cout<<"Polinomul nu are radacini in intervalul dat.";
            clearscreen();
        }

        else if (choice=='3')
        {
            Polinom<NrRational> f,g;
            alg<NrRational>(f,g);
            clearscreen();
        }
        else if (choice=='4')
        {
            Polinom<NrComplex> f,g;
            alg<NrComplex>(f,g);
            clearscreen();
        }
        else if (choice=='x' || choice=='X')
            break;
        else
        {
            cout<<"Varianta nu exista. Apasati X pentru a incerca din nou.\n";
            char ws;
            cin>>ws;
            if (ws=='x' || ws=='X')
                system("CLS");
        }
    }
    return 0;
}

/*
int
3
2 3 0 -1
2
3 2 1

float
2
1 -2 1
3
2.1 3 0 -1.4

NrRational
2
2 1 3 2 1 1
1
3 1 2 4

NrComplex
2
2 1 3 2 1 1
1
3 1 2 4
*/
