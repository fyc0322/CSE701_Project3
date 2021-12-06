#include <string>
#include <stack>
#include <cmath>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
using namespace std;
const string::size_type ptr_num=9;  // number of operators.
const string::size_type fun_num=24;  // number of functions.
string Opr_type[] = {"Limit1", "Integration", "Fourier_coefficients", "Get_zero", "Integration_0_inf", "FFT", "random","value"}; //calculations can be done.
// types of words in the input expression.
typedef enum
{
    TKT_NUMBER,    // number 
    TKT_OPERATOR,  // operator 
    TKT_FUNCTION,  // function 
    TKT_ENDSIGN,   // Terminator 
    TKT_UNKNOW     // Unknown symbol 
} TokenType;


class Calculator
{
public:
    Calculator(const string &expr,const string &deri); // initializing Calculator with a string.
    Calculator(const Calculator &value); // copy constructor.
    void print();   // print expr
    void update(); // update the length of expr.
    bool num_detect(string str); // detect the type of input.
    void getVal(double &res); // res is the solution value.
    double cal_num(const string &cal_expr, double x); // get an arithmetic value of expr.
    double integra(string& expr_i, double a, double b, double eps); // do integration.
    double integra_0i(string &expr_0i); // do integration from 0 to infinity.
    vector<double> zero(string &expr0,double a, double b, double h, double eps); // get zero value.
    void limit1(double a, double b, double h, double eps); // get local limit value with f(x) and f'(x).
    void fourier(vector<double>::size_type n); // get fourier coefficients
    void fft(vector<double>::size_type k,vector<double>::size_type il,double h); // fast fourier transform.
    void random(double r,double u, double g, vector<double>::size_type n);
    void value(double a, double b,double h);


private:
    string expr;   // store the mathematical expression for calculation.
    string orig;   // store original mathematical expression.
    string deri;  // store the derivative expression of orig.
    string token; // the word read at each time. 
    TokenType tkType;  // word type. 
    string::size_type pos, length;   // position for read and its length. 

    static string ptrList[];   // list of supported operators.
    static int32_t ptrArgCnt[];         // the number of arguments required for the operator.
    static string funList[];   // list of supported functions. 
    static int32_t funArgCnt[];         // the number of arguments required for the function. 
    static int32_t preceMap[][ptr_num]; // operation priority table.

    // read next word
    void readToken();

    // compare the priority of the two operators (ptr1 and ptr2), and then return an int32_t value.
    int32_t comparePrece(const string &ptr1, const string &ptr2);

    // single-step calculation of an operator and return a double value.
    double calculate(const string &ptr, double arg[]);

    // single-step calculation of a function and return a double value.
    double callFun(const string &fun, double arg[]);

    // get operator (ptr) sequence number, and then return an int32_t value.
    int32_t getPtrIndex(const string &ptr);

    // get function (fun) sequence number, and then return an int32_t value.
    int32_t getFunIndex(const string &fun);

    // check whether the number of arguments in opnd matches n.
    // if match, get n parameters from operand stack (opnd) and store them in arg.
    bool getArg(stack<double> &opnd, double arg[], int32_t n);
};

// initialize static data. "f(" means function and "#" is the terminator.
string Calculator::ptrList[] = {"f(", ",", "+", "-", "*", "/", "(", ")", "#"};

int32_t Calculator::ptrArgCnt[] = {0, 1, 2, 2, 2, 2, 0, 0, 0};

string Calculator::funList[] = 
{
    "sin(",  "cos(",  "tan(",
    "acos(", "asin(", "atan(",
    "cosh(", "sinh(", "tanh(",
    "acosh(","asinh(","atanh(",
    "log(",  "log2(", "log10(", 
    "exp(",  "exp2(", "sqrt(",
    "pow(",  "cbrt(", "hypot(",
    "fmax(", "fmin(", "abs(",
};

int32_t Calculator::funArgCnt[] = 
{
    1, 1, 1,
    1, 1, 1,
    1, 1, 1,
    1, 1, 1,
    1, 1, 1,
    1, 1, 1,
    2, 1, 2,
    2, 2, 1,
};

// defination of preceMap is illustrated in external doc.
int32_t Calculator::preceMap[][ptr_num] = 
{
  //{f(, , +, -, *, /, (, ), #}
    {-1,-1,-1,-1,-1,-1,-1, 1, 1},  // f(
    { 1, 1, 1, 1, 1, 1, 1, 1, 1},  // ,
    {-1, 1, 1, 1,-1,-1,-1, 1, 1},  // +
    {-1, 1, 1, 1,-1,-1,-1, 1, 1},  // -
    {-1, 1, 1, 1, 1, 1,-1, 1, 1},  // *
    {-1, 1, 1, 1, 1, 1,-1, 1, 1},  // /
    {-1, 1,-1,-1,-1,-1,-1, 0, 1},  // (
    { 1, 1, 1, 1, 1, 1, 1, 1, 1},  // )
    {-1, 1,-1,-1,-1,-1,-1,-1, 0},  // #
};

/******************************************************************************/

Calculator::Calculator(const string &expr,const string &deri)
{
    this->expr = expr + "#";
    this->orig = expr + "#";
    this->deri = deri + "#";
    this->token = "";
    this->pos = 0;
    this->length = 1 + expr.length();
}

Calculator::Calculator(const Calculator &value)
{
    this->expr = value.expr;
    this->orig = value.orig;
    this->deri = value.deri;
    this->token = value.token;
    this->pos = 0;
    this->length = value.length;
}

void Calculator::print()
{
    cout << orig << endl;
    cout << deri << endl;
}

void Calculator::update()
{
    length = expr.length();
}

//detect whether str is a number of a string of letters.
bool Calculator::num_detect(string fstr)
{
    stringstream str_(fstr);
    double d;
    string c;
    if(!(str_ >> d))
        return false;
    if (str_ >> c)
        return false;
    return true;
}

int32_t Calculator::getPtrIndex(const string &ptr)
{
    for (string::size_type i = 0; i < ptr_num; ++ i)
    {
        if (ptrList[i] == ptr) return static_cast<int32_t>(i);
    }
    return -1;
}

int32_t Calculator::getFunIndex(const string &fun)
{
    for (string::size_type i = 0; i < fun_num; ++ i) 
    {
        if (funList[i] == fun) return static_cast<int32_t>(i);
    }
    return -1;
}

int32_t Calculator::comparePrece(const string &ptr1, const string &ptr2)
{
    int32_t m = getPtrIndex(ptr1);
    int32_t n = getPtrIndex(ptr2);
// if m or n == -1, then m or n is a function.
    if (m == -1) m = 0;
    if (n == -1)n = 0;
    return preceMap[m][n];
}

double Calculator::calculate(const string &ptr, double arg[])
{
    switch (getPtrIndex(ptr)) 
    {
        case 1:
            return arg[0];

        case 2:
            return arg[0] + arg[1];

        case 3:
            return arg[0] - arg[1];

        case 4:
            return arg[0] * arg[1];

        case 5:
            return arg[0] / arg[1];
    }
    return 0;
}

double Calculator::callFun(const string &fun, double arg[])
{
    switch(getFunIndex(fun)) 
    {
        case 0:
            return sin(arg[0]);

        case 1:
            return cos(arg[0]);

        case 2:
            return tan(arg[0]);

        case 3:
            if (arg[0] > 1 || arg[0]<-1)
            {
                cout << "Error: invalid input for acos().\n" << endl;
                exit(3);
            }
            return acos(arg[0]);

        case 4:
            if (arg[0] > 1 || arg[0]<-1)
            {
                cout << "Error: invalid input for asin().\n" << endl;
                exit(3);
            }
            return asin(arg[0]);

        case 5:
            return atan(arg[0]);

        case 6:
            if (cosh(arg[0]) > 1e15)
            {
                cout << "Error: invalid input for cosh().\n" << endl;
                exit(3);
            }
            return cosh(arg[0]);

        case 7:
            if (sinh(arg[0]) > 1e15)
            {
                cout << "Error: invalid input for sinh().\n" << endl;
                exit(3);
            }
            return sinh(arg[0]);

        case 8:
            if (tanh(arg[0]) < 1e-15)
            {
                cout << "Error: invalid input for tanh().\n" << endl;
                exit(3);
            }
            return tanh(arg[0]);

        case 9:
            if (arg[0] < 1)
            {
                cout << "Error: invalid input for acosh().\n" << endl;
                exit(3);
            }
            return acosh(arg[0]);

        case 10:
            return asinh(arg[0]);

        case 11:
            if (arg[0] >= 1 || arg[0]<=-1)
            {
                cout << "Error: invalid input for atanh().\n" << endl;
                exit(3);
            }
            return atanh(arg[0]);

        case 12:
            if ( arg[0] <= 0)
            {
                cout << "Error: invalid input for log().\n" << endl;
                exit(3);
            }
            return log(arg[0]);

        case 13:
            if ( arg[0] <= 0)
            {
                cout << "Error: invalid input for log2().\n" << endl;
                exit(3);
            }
            return log2(arg[0]);

        case 14:
            if ( arg[0] <= 0)
            {
                cout << "Error: invalid input for log10().\n" << endl;
                exit(3);
            }
            return log10(arg[0]);

        case 15:
            if ( exp(arg[0]) > 1e15 || exp(arg[0]) < 1e-15)
            {
                cout << "Error: invalid input for exp().\n" << endl;
                exit(3);
            }
            return exp(arg[0]);

        case 16:
            if ( exp2(arg[0]) > 1e15 || exp2(arg[0]) < 1e-15)
            {
                cout << "Error: invalid input for exp2().\n" << endl;
                exit(3);
            }
            return exp2(arg[0]);

        case 17:
            if ( arg[0] < 0)
            {
                cout << "Error: invalid input for sqrt().\n" << endl;
                exit(3);
            }
            return sqrt(arg[0]);

        case 18:
            return pow(arg[0],arg[1]);

        case 19:
            return cbrt(arg[0]);

        case 20:
            if ( arg[0] <= 0 || arg[1] <= 0)
            {
                cout << "Error: invalid input for hypot().\n" << endl;
                exit(3);
            }
            return hypot(arg[0],arg[1]);

        case 21:
            return fmax(arg[0],arg[1]);

        case 22:
            return fmin(arg[0],arg[1]);

        case 23:
            return abs(arg[0]);

    }
    return 0;
}

void Calculator::readToken()
{
// check whether the reading is completed.
    if (pos >= length) 
    {
        tkType = TKT_ENDSIGN;
        return ;
    }
// record the current read position and corresponding characters
    string::size_type pos_t = pos;
    char ch = expr[pos_t++];
// judge whether ch is an operator or negative number
    if (-1 != getPtrIndex(std::string(1, ch))) 
    {
        if (ch != '-') 
        {
            tkType = TKT_OPERATOR;
        }
        else 
        {
            if (pos_t>1 && (expr[pos_t-2] == ')' || isdigit(expr[pos_t-2]))) 
            {
                tkType = TKT_OPERATOR;
            }
            else
            {
                ++pos_t;
                while (pos_t < length && isdigit(ch = expr[pos_t])) 
                {
                    ++pos_t;
                }

                if (ch == '.') 
                {
                    ++pos_t;
                    while (pos_t < length && isdigit(expr[pos_t]))
                {
                    ++pos_t;
                }
                }
                tkType = TKT_NUMBER;
            }
        }
    }
// judge whether ch is a number
    else if (isdigit(ch)) 
    {
        while (pos_t < length && isdigit(ch = expr[pos_t])) 
        {
            ++pos_t;
        }

        if (ch == '.') 
        {
            ++pos_t;
            while (pos_t < length && isdigit(expr[pos_t])) ++pos_t;
        }
        tkType = TKT_NUMBER;
    } 
// judge whether ch is a function
    else if (isalpha(ch)) 
    {
        while (pos_t < length && (isalnum(ch) || ch == '_')) 
        {
            ch = expr[++pos_t];
        }

        if (ch == '(') 
        {
            ++pos_t;
            if (-1 != getFunIndex(expr.substr(pos, pos_t - pos))) tkType = TKT_FUNCTION;
            else tkType = TKT_UNKNOW;
        }
        else tkType = TKT_UNKNOW;
    } 
    else tkType = TKT_UNKNOW;
    token = expr.substr(pos, pos_t - pos);
 //   cout << token <<endl;
    pos = pos_t;
}

bool Calculator::getArg(stack<double> &opnd, double arg[], int32_t n)
{
    if (opnd.size() < static_cast<unsigned long>(n)) 
    {
        return false;
    }

    for (int32_t i = n - 1; i >= 0; -- i) 
    {
        arg[i] = opnd.top();
        opnd.pop();
    }
    return true;
}

void Calculator::getVal(double &res)
{
    stack<string> optr;  // stack for operators
    stack<double> opnd;       // stack for operands 
    int32_t comRes;                    // store priority of two operators.
    int32_t idx;                       // store sequence number of a function or an operator.
    int32_t argCnt;                    // store the required number to operate a function or an operator.
    optr.push("#");
    pos = 0;
    readToken();
// calculate the value of input expression.
    while(tkType != TKT_ENDSIGN || !optr.empty()) 
    {

        if(tkType == TKT_UNKNOW)
        {
            cout << "Error: The input mathematical expression contains unknown characters\n" << endl;
            exit(1);
        }
// if the word from readToken is a number, then push it in opnd.
        if (tkType == TKT_NUMBER) 
        {
            opnd.push(atof(token.c_str()));
            readToken();
        }
// do the calculation based on operation priority.
        else 
        {
            comRes = comparePrece(optr.top(), token);
            switch(comRes) 
            {
                case -1:
                    optr.push(token);
                    readToken();
                    break;

                case 1:
                {
                    string ptr = optr.top();
                    optr.pop();
                    idx = getPtrIndex(ptr);
                    double arg[ptr_num]; 
                    if (-1 != idx) 
                    {
                        argCnt = ptrArgCnt[idx];
                        if (argCnt) 
                        {
                            if (!getArg(opnd, arg, argCnt))
                            {
                                cout << "Error: The number of operands failed to match the requirement\n" << endl;
                                exit(2);
                            }
                            res = calculate(ptr, arg);
                            opnd.push(res);
                        }
                    } 
                    else 
                    {
                        idx = getFunIndex(ptr);
                        argCnt = funArgCnt[idx];
                        if (!getArg(opnd, arg, argCnt))
                        {
                            cout << "Error: The number of operands failed to match the requirement\n" << endl;
                            exit(2);
                        }
                        res = callFun(ptr, arg);
                        opnd.push(res);
                        readToken();
                    }
                    break;
                }

                case 0:
                    optr.pop();
                    readToken();
                    break;
            }
        }
    }
    res = opnd.top();
}

double Calculator::cal_num(const string &cal_expr,double x)
{
    double res;
    if (abs(x)<1e-8) 
    {
        x=0;
    }
    stringstream ss;
    if (cal_expr.find('#')== string::npos) 
    {
        expr=cal_expr + "#";
    }
    else expr=cal_expr;
    ss << setprecision(5) << x;
    string str_n=ss.str();
    while (expr.find('X') < expr.length()) 
    {
        expr.replace(expr.find('X'), 1, str_n);
    }
    update();
    getVal(res);
    return res;
}

double Calculator::integra(string& expr_i, double a, double b, double eps){
// n: divide the integral interval into 2n equal parts 
// fa and fb are the lower and upper boundary.
// h is the width of the integral interval
// t1 is the return value of integration result.
// other variables are for temporary usage.
    int32_t n,k;
    double fa,fb,h,t1,p,s,x,t=10,l;
    if (a>b)
    {
        cout<<"Integration Error: upper limit is lower than lower limit: " <<endl;
        exit(2);
    }

// get function value at left and right boundary
    fa=cal_num(expr_i,a);
    fb=cal_num(expr_i,b);
// do the integration
    n=1;
    h=b-a;
    t1=h*(fa+fb)/2.0;
    p=1.0;
    while (p>=eps)
    {
        s=0.0;
        for (k=0;k<n;k++)
        {
            x=a+(k+0.5)*h;
            l=cal_num(expr_i,x);
            s=s+l;
        }
        t=(t1+h*s)/2.0;
        p=fabs(t1-t);
        t1=t; 
        n=n+n; 
        h=h/2.0;
    }
    return(t);
  }

double Calculator::integra_0i(string &expr_0i)
{
// x are the roots of Laguerre polynomial Ln(x)
// lc are the weights for integration
// g is the integration result 
    vector<double>::size_type i;
    double x,g;
    vector<double> t{0.26355990,1.41340290,3.59642600,7.08580990,12.64080000};
    vector<double> Lc{0.6790941054,1.638487956,2.769426772,4.315944000,7.104896230};
    g=0.0;
    for (i=0; i<=4; i++)
    { 
        x=t[i]; 
        g=g+Lc[i]*cal_num(expr_0i,x); 
    }
    return(g);

}

vector<double> Calculator::zero(string &expr0,double a, double b, double h, double eps)
{
// a and b are interval endpoints.
// h is step length and eps is the precision.
// n is the number of zero-points.
// x is a group of zero points.
// y, y0 and y1 are the function values at end points.
// z, z0 and z1 are  end points.
    int n,js;
    vector<double> x;
    double z,y,z1,y1,z0,y0;
    n=0; 
    z=a; 
    y=cal_num(expr0,z);
    while ((z<=b+h/2.0))
    { 
        if (fabs(y)<eps)
        { 
            n=n+1; 
            x.push_back(z);
            z=z+h/2.0; 
            y=cal_num(expr0,z);
        }
        else
        { 
            z1=z+h; 
            y1=cal_num(expr0,z1);
            if (fabs(y1)<eps)
            {
                n=n+1; 
                x.push_back(z1);
                z=z1+h/2.0; 
                y=cal_num(expr0,z);
            }
            else if (y*y1>0.0)
            {
                y=y1; 
                z=z1;
            }
            else
            { 
                js=0;
                while (js==0)
                {
                    if (fabs(z1-z)<eps)
                    {
                        n=n+1; 
                        x.push_back((z1+z)/2.0);
                        z=z1+h/2.0; 
                        y=cal_num(expr0,z);
                        js=1;
                    }
                    else
                    {
                        z0=(z1+z)/2.0;
                        y0=cal_num(expr0,z0);
                        if (fabs(y0)<eps)
                        {
                            x.push_back(z0); 
                            n=n+1; 
                            js=1;
                            z=z0+h/2.0; 
                            y=cal_num(expr0,z);
                        }
                        else if ((y*y0)<0.0)
                        {
                            z1=z0; 
                            y1=y0;
                        }
                        else 
                        {
                            z=z0; 
                            y=y0;
                        }
                    }
                }
            }
        }
    }
    for (vector<double>::size_type i=0; i<x.size(); i++)
        if (abs(x[i])<1e-10) 
        {
            x[i]=0;
        }
    if (n==0) 
    {
        cout << "function at interval " << a << " to "<< b << " has no zero value.\n" << endl;
    }
    if (n>0){
        cout << "x-points for 0 from " << a << " to "<< b << " are: ";
        for(vector<double>::size_type i=0; i<x.size();i++) 
        {
            cout << x[i] << " " ;
        }
        cout << "\n" <<endl;
    }
    return(x);
}

void Calculator::limit1(double a, double b, double h, double eps){
    vector<double> dr_zero;
    vector<double>::size_type n;
    double lb,rb;
    dr_zero=zero(deri,a,b,h,eps);
    n=dr_zero.size();
    if (n == 0) 
    {
        cout << "no limit value between " << a << " to "<< b << ". largest and lowest value are both on the boundaries"<<endl;
        lb=cal_num(orig,a);
        rb=cal_num(orig,b);
        (lb>rb) ? cout << "highest value is: " << lb << " lowest value is: " << rb <<"\n"<<endl :
                  cout << "highest value is: " << rb << " lowest value is: " << lb <<"\n"<<endl;

    }
    else{
        for(vector<double>::size_type i=0; i<dr_zero.size();i++)
        {
            lb = cal_num(orig,dr_zero[i]);
            rb = cal_num(orig,dr_zero[i]+0.1);
            (lb > rb) ? cout << dr_zero[i] << " is the local maxima point. The maxima is: " << lb <<endl :
                        cout << dr_zero[i] << " is the local minimum point. The minimum is: " << lb <<"\n"<<endl;
        }
    }
}


void Calculator::fourier(vector<double>::size_type n)
{
// a is for Fourier's coefficient ak.
// b is for Fourier's coefficient bk.
// f is the function value at 2n+1 points.
// t is the constant.
// c and s are cos(2*pi/(2n+1)) and sin(2*pi/(2n+1))
// c1 and s1 are cos(kθ) and sin(kθ)
// u1,u2,u0 are temporary variables.
    vector<double>::size_type i,j;
    vector<double> a,b,f;
    double t,c,s,c1,s1,u1,u2,u0;
    t=6.283185306/(2.0*n+1.0);
    c=cos(t); 
    s=sin(t);
    for (i=0; i<=2*n; i++) 
    {
        f.push_back(cal_num(orig,(i+0.5)*t));
    }
    t=2.0/(2.0*n+1.0); 
    c1=1.0; 
    s1=0.0;
    for (i=0; i<=n; i++)
    { 
        u1=0.0; 
        u2=0.0;
        for (j=2*n; j>=1; j--)
        { 
            u0=f[j]+2.0*c1*u1-u2;
            u2=u1; 
            u1=u0;
        }
        a.push_back(t*(f[0]+u1*c1-u2));
        b.push_back(t*u1*s1);
        u0=c*c1-s*s1; 
        s1=c*s1+s*c1; 
        c1=u0;
    }
// output operation
    ofstream output;
    output.open ("fourier_coefficients.txt");
    if (!output)
    {
        cout<<"\nError: failed to open the input file." <<endl;
        exit(0);
    }
    else
    {
        output << "Math formula is: " << orig << "\n"<< endl;
        output << " i         " <<"a           " <<"   b"<<endl;
        for(vector<double>::size_type i=0; i<a.size();i++)
        {
            output <<setw(2)<< i << "      "<< setw(11)<<a[i] <<"       "<< b[i] << endl;
        }
    }
    output.close();
    return;
}

void Calculator::fft(vector<double>::size_type k,vector<double>::size_type il,double h)
{
// pr is used to store the real part of input.
// pi is used to store the imaginary part of input.
// fr is used to return the real part of result.
// fi is used to return the imaginary part of result.
// n is the sample points.
// h is the sample starting points.
// il = 0 means not to calculate magnitude and angle; = 1 means to calculate them.
// it,m,is,i,j,nv,l0 are all for loop-iteration.
// p,q,s,tr,ti,pdr,pdi are all for temporary purpose.
    vector<double>::size_type n = static_cast<vector<double>::size_type>(pow(2,k));
    vector<double> pr(n,0),pi(n,0),fr(n,0),fi(n,0),sq;
    vector<double>::size_type it,m,is,i,j,nv,l0; // for iteration
    double p,q,s,tr,ti,pdr,pdi;

// get value of original points.
    for (i=0; i<n; i++)
        {
            pr[i]=cal_num(orig,0.1*i+0.1*h); 
            pi[i]=0.0;
        }
    sq=pr;

// get corresponding Even power and odd power.
    for (it=0; it<n; it++)
        { 
            m=it; 
            is=0;
            for (i=0; i<k; i++)
            { 
                j=m/2; 
                is=2*is+(m-2*j); 
                m=j;
            }
            fr[it]=pr[is]; 
            fi[it]=pi[is];
        }

// do FFT to get the results
        pr[0]=1.0; 
        pi[0]=0.0;
        p=6.283185306/(1.0*n);
        pr[1]=cos(p); 
        pi[1]=-sin(p);
        for (i=2; i<n; i++)
        { 
            p=pr[i-1]*pr[1]; 
            q=pi[i-1]*pi[1];
            s=(pr[i-1]+pi[i-1])*(pr[1]+pi[1]);
            pr[i]=p-q; 
            pi[i]=s-p-q;
        }
        for (it=0; it<n-1; it=it+2)
        { 
            tr=fr[it]; 
            ti=fi[it];
            fr[it]=tr+fr[it+1]; 
            fi[it]=ti+fi[it+1];
            fr[it+1]=tr-fr[it+1]; 
            fi[it+1]=ti-fi[it+1];
        }
        m=n/2; 
        nv=2;
        for (l0=k-1; l0>=1; l0--)
        { 
            m=m/2; 
            nv=2*nv;
            for (it=0; it<=(m-1)*nv; it=it+nv)
                for (j=0; j<=(nv/2)-1; j++)
                { 
                    p=pr[m*j]*fr[it+j+nv/2];
                    q=pi[m*j]*fi[it+j+nv/2];
                    s=pr[m*j]+pi[m*j];
                    s=s*(fr[it+j+nv/2]+fi[it+j+nv/2]);
                    pdr=p-q; 
                    pdi=s-p-q;
                    fr[it+j+nv/2]=fr[it+j]-pdr;
                    fi[it+j+nv/2]=fi[it+j]-pdi;
                    fr[it+j]=fr[it+j]+pdr;
                    fi[it+j]=fi[it+j]+pdi;
                }
        }
//output operation
    if (il==1)
        for (i=0; i<n; i++)
        { 
            pr[i]=sqrt(fr[i]*fr[i]+fi[i]*fi[i]);
            if (fabs(fr[i])<0.000001*fabs(fi[i]))
            { 
                if ((fi[i]*fr[i])>0) pi[i]=90.0;
                else pi[i]=-90.0;
            }
            else
                pi[i]=atan(fi[i]/fr[i])*360.0/6.283185306;
        }

    ofstream output;
    output.open ("FFT_result.txt");
    if (!output)
    {
        cout<<"\nError: failed to open the input file." <<endl;
        exit(0);
    }
    else
    {
        output << "Math formula is: " << orig << "\n"<< endl;
        output << " f(X)       " << setw(15) <<"\"Real part after FFT\"" <<"   Imag part after FFT "<<endl;
        for(vector<double>::size_type i=0; i<n;i++)
        {
            output <<setw(10)<< sq[i] << setw(15)<< fr[i] <<"       "<< setw(15)<<fi[i] << endl;
        }
        output<<" "<<endl;

        if (il == 1)
            output << " f(X)       " << setw(15) <<"\"magnitude\"" <<"               Angle"<<endl;
            for(vector<double>::size_type i=0; i<n;i++)
            {
            output <<setw(10)<< sq[i] << setw(15)<< pr[i] <<"       "<< setw(15)<<pi[i] << endl;
            }
    }
    output.close();
    return;
}

void Calculator::random(double r,double u, double g, vector<double>::size_type n)
{
// r is the random seed.
// u is the average value.
// g is the variance of normal distribution. when g<=0, we ignore normal distribution restriction.
// n is the length of data
// s, w, v are function coefficients listed in doc.
// t and m for temporary usage.
    vector<double>::size_type i,k;
    vector<double> data(n,0);
    int m;
    double s,w,v,t;
    s=65536.0; 
    w=2053.0; 
    v=13849.0;
    if (g > 0)
    {
        for (k=0; k<n; k++)
        {
            t=0.0;
            for (i=1; i<=12; i++)
            { 
                r=r*w+v;
                m=(int)(r/s);
                r=r-m*s;
                t=t+r/s;
            }
            data[k]=u+g*(t-6.0);
        }
    }
    else
    {
        for (k=0; k<n; k++)
        {
            m=(int)(r/s); 
            r=r-m*s;
            r=w*r+v; 
            m=(int)(r/s);
            r=r-m*s; 
            data[k]=r/s;
        }
    }
// output operation
    ofstream output;
    output.open ("random.txt");

    if (!output)
    {
        cout<<"\nError: failed to open the input file." <<endl;
        exit(0);
    }

    else
    {
        if (g > 0)
        {
            output << "Normal distribution random points with average= "<< u << " variance= "<< g <<":"<<"\n"<< endl;
            for(vector<double>::size_type i=0; i<data.size();i++)
            {
                output <<setw(2)<< data[i] << ",  ";
            }
            output << "End" << endl;
        }
        else
        {
            output << "random points in [0,1] are: " <<"\n"<< endl;
            for(vector<double>::size_type i=0; i<data.size();i++)
            {
                output << data[i] << ",  ";
            }
            output << "End" << endl;
        }
    }
    output.close();
    return;
}

void Calculator::value(double a, double b,double h)
{
    double i=a;
    vector<double> t,ti;
    while(i < b)
    {
        t.push_back(cal_num(orig,i));
        ti.push_back(i);
        i=i+h;
    }
    for (vector<double>::size_type j=0;j<t.size();j++)
    {
        cout << "f(" << ti[j] <<")= " << t[j]<<endl;
    }
}

void cal_expression()
{
    string t;         // store temporary data
    string str;       // store a whole row from input.txt
    string Math_expr; // store mathematical expression
    string Deri_expr = "none"; // store derivative expression.
    string str_tem;   // store the operation of Math_expr 
    double res;       // arithmetic calculation result
    vector<double> parameter,zeros;
    int32_t type=99;  // identify which operation in Opr_type is going to operate.
    int32_t Opr_num = end(Opr_type) - begin(Opr_type); // length of Opr_type
// read mathematical expression from input.txt
    ifstream input("input.txt");
// detect whether the file is open or not.
    if(!input.is_open())
    {
    cout << "\nError: failed to open the input file." << endl;
    exit(1);
    }

// first loop to read the expression.
    while(getline(input,str))
    {

// remove spaces in str
        str.erase(remove(str.begin(), str.end(), ' '), str.end());

// extract original expression from file into Math_expre
        if(str.find("f(X)")<str.length() && str.find('#')== string::npos)
        {
            Math_expr = str.substr(str.find("f(X)")+5);
            cout << "f(X)=" << Math_expr <<endl;

// if there is no variable in expression, then just do arithmetic calculation
            if (Math_expr.find('X')== string::npos)
            {
                Calculator cal(Math_expr,Deri_expr);
                cal.getVal(res);
                cout << Math_expr << "=" << res << endl;
                return;
            }
        }
// extract derivative expression of expr
        if (str.find("f'(X)")<str.length() && str.find('#')== string::npos)
        {
            Deri_expr = str.substr(str.find("f'(X)")+6);
            cout << "f'(X)=" << Deri_expr <<"\n"<<endl;
        }
        if(str.find("Arithmatic")<str.length()) break;
    }

    Calculator cal(Math_expr,Deri_expr);

// second loop to do specific mathematical operation.
    while(getline(input,str))
    {
        type=99;
        if(str.find("do:")<str.length() && str.find('#')== string::npos) 
        {
            str_tem = str.substr(str.find("do:")+3);
            stringstream tem(str_tem);
            tem >> t;
            for (int32_t i = 0; i < Opr_num; i++) 
            {
                if (Opr_type[i] == t) type = i;
            }

            switch (type) 
            {
                case 0:
                    cout << "***** Limit Operation *****"<<endl;
                    while(tem >>t)
                    {
                        if (!cal.num_detect(t))
                        {
                            cout << "\nError: invalid input format for limit. It should be a number."<<endl;
                            exit(2);
                        }
                        parameter.push_back(stod(t));
                    }
                    if (parameter.size() != 4)
                    {
                        cout <<"*****Error: too many or not enough data for limit.*****" << endl;
                        exit(2);
                    } 
                    cal.limit1(parameter[0],parameter[1],parameter[2],parameter[3]);
                    parameter.clear();
                    break;

                case 1:
                    cout << "***** Integration Operation *****"<<endl;
                    // vector<double> parameter [0] and [1] are lower and upper limit of integration, [2] is the precision. 
                    while(tem >>t)
                    {
                        if (!cal.num_detect(t))
                        {
                            cout << "\nError: invalid input format for integration. It should be a number."<<endl;
                            exit(2);
                        }
                        parameter.push_back(stod(t));
                    }
                    if (parameter.size() != 3)
                    {
                        cout <<"Error: too many or not enough data for integration." << endl;
                        exit(2);
                    } 
                    parameter.push_back(cal.integra(Math_expr,parameter[0],parameter[1],parameter[2]));
                    cout<<"integration from "<<parameter[0]<< " to "<< parameter[1] << " of " << Math_expr << " is: "<<parameter[3]<<"\n"<<endl;
                    parameter.clear();
                    break;

                case 2:
                    cout << "***** Fourier_coefficients Operation *****" << endl;
                    while(tem >>t)
                    {
                        if (!cal.num_detect(t))
                        {
                            cout << "\nError: invalid input format for Fourier_transform. It should be a number."<<endl;
                            exit(2);
                        }
                        parameter.push_back(stod(t));
                    }
                    if (parameter.size() != 1)
                    {
                        cout <<"Error: too many or not enough data for Fourier_transform." << endl;
                        exit(2);
                    }
                    if (parameter[0] < 0)
                    {
                        cout <<"Error: input number for Fourier_transform should be a positive number." << endl;
                        exit(2);
                    }
                    cal.fourier(static_cast<vector<double>::size_type>(parameter[0]));
                    cout << "Successfully calculate Fourier_coefficients and store results in fourier_coefficients.txt.\n";
                    parameter.clear();
                    break;

                case 3:
                    cout << "***** Get zero Operation *****"<<endl;
                    while(tem >>t)
                    {
                        if (!cal.num_detect(t))
                        {
                            cout << "\nError: invalid input format for Get zero. It should be a number."<<endl;
                            exit(2);
                        }
                        parameter.push_back(stod(t));
                    }
                    if (parameter.size() != 4)
                    {
                        cout <<"Error: too many or not enough data for Get zero." << endl;
                        exit(2);
                    } 
                    zeros=cal.zero(Math_expr,parameter[0],parameter[1],parameter[2],parameter[3]);
                    parameter.clear();
                    break;

                case 4:
                    cout << "***** Integrate from 0 to infinity *****"<<endl;
                    while(tem >>t) 
                    {
                        cout <<"Error: too many arguments for integration from 0 to infinity." << endl;
                        exit(2);
                    }
                    cout << "integration from 0" << " to "<< " infinity of " << Math_expr << " is: "<<cal.integra_0i(Math_expr)<<"\n"<<endl;
                    break;

                case 5:
                    cout << "\n***** FFT Operation *****"<<endl;
                    while(tem >>t)
                    {
                        if (!cal.num_detect(t))
                        {
                            cout << "\nError: invalid input format for Get zero. It should be a number."<<endl;
                            exit(2);
                        }
                        parameter.push_back(stod(t));
                    }
                    if (parameter.size() != 3)
                    {
                        cout <<"Error: too many or not enough data for Get zero." << endl;
                        exit(2);
                    }
                    if (parameter[0]<=1 || !(parameter[1] == 0 || parameter[1] == 1) || parameter[2] == 0)
                    {
                        cout << "error: input for FFT have not suitable value." << endl;
                        exit(3);
                    }
                    cal.fft(static_cast<vector<double>::size_type>(parameter[0]),static_cast<vector<double>::size_type>(parameter[1]),parameter[2]);
                    cout << "Successfully calculate FFT and store results in FFT_result.txt.\n";
                    parameter.clear();
                    break;

                case 6:
                    cout << "\n***** random Operation *****"<<endl;
                    while(tem >>t)
                    {
                        if (!cal.num_detect(t))
                        {
                            cout << "\nError: invalid input format for random. It should be a number."<<endl;
                            exit(2);
                        }
                        parameter.push_back(stod(t));
                    }

                    if (parameter.size() != 4)
                    {
                        cout <<"Error: too many or not enough data for Get zero." << endl;
                        exit(2);
                    }

                    if (parameter[3]<1)
                    {
                        cout << "error: input for number of random data should be positive." << endl;
                        exit(3);
                    }
                    cal.random(parameter[0],parameter[1],parameter[2],static_cast<vector<double>::size_type>(parameter[3]));
                    cout << "Successfully calculate random points and store results in random.txt.\n";
                    parameter.clear();
                    break;

                case 7:
                    cout << "\n***** value Operation *****"<<endl;
                    while(tem >>t)
                    {
                        if (!cal.num_detect(t))
                        {
                            cout << "\nError: invalid input format for calculating value. It should be a number."<<endl;
                            exit(2);
                        }
                        parameter.push_back(stod(t));
                    }

                    if (parameter.size() != 3)
                    {
                        cout <<"Error: too many or not enough data for value." << endl;
                        exit(2);
                    }

                    if (parameter[0]>parameter[1] || parameter[2]<=0)
                    {
                        cout << "error: a should be less than b and st should be positive." << endl;
                        exit(3);
                    }
                    cal.value(parameter[0],parameter[1],parameter[2]);
                    parameter.clear();
                    break;

                case 99:
                    cout << "Error: invalid input to select mathematical operation"<<endl;
                    exit(2);
            } 
        }
    }
    input.close();
}








