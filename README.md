# CSE701_Project3 - Mathematic Calculator

## Author name: 

Yucheng Feng

## Mission:

Read mathematical expressions from text file. Based on the formula loaded, the program should achieve goals of performing popular arithmetic calculations, such as + - * / common functions in “cmath”, doing integration, calculating Fourier coefficients, getting zero values of expressions, finding limit values of the expression in the given range, doing FFT and creating random numbers.

## Introduction:
“Mathematic Calculator” is fulfilled with C++ language and designed for mathematical operation in production-level. It can read formula from “input.txt” and then do sin, cos, tan, acos, asin, atan, cosh, sinh, tanh, acosh, asinh, atanh, log, log2, log10, exp, exp2, sqrt, pow, cbrt, hypot, fmax, fmin, abs, which are quite popular functions in C++ library “cmath”. It can also calculate integration, find limit-value, get zero-value, get Fourier coefficients, do FFT and create random points.

Codes for “Mathematic Calculator” are written in two files – calculator.h and Math_co.cpp. Most jobs are done in calculator.h, in which corresponding functions are encapsulated in Class “Calculator”. Math_co.cpp is just used to call a function to read formulas, print out or save calculation results.

“Mathematic Calculator” works for arbitrary of formulas, but it has its own limitations which will be described in “Cautions” part of this document. Users of this code just need to give their formula in input.txt file, then run Math_co.cpp to get the results they are looking forward to.

## Input:
All input information should be listed in “input.txt”. This file has two parts. The first part is designed to read the formula and the second is used to do specific calculation. To read formulas from input file, users should give their own mathematical expressions after “The input mathematical expression is:” and “The input mathematical expression's derivative:”. These sentences cannot be changed otherwise it will lead to read file errors. The format of basic functions should be in the same format as those in the C++ cmath library. The second is to select specific operations. Nowadays, this code offers eight possible operations; users should choose their own operations based on their needs, and typed them after “The operation you wanna do:”. This sentence cannot be changed again.The format of each operation is listed below.

**Integration a b eps**: where Integration indicates that we are going to do integration from a to b with precision eps.

**Get_zero a b h eps**: where Get_zero indicates that we are going to find all the zero values of input formula in the range of a and b with step h and precision eps. Get_zero is based on Bisection method so it needs a step.

**Limit1 a b  h  eps**: where Limit1 indicates that we are going to find the limit value of formula in given range and other input parameters’ meaning are just like these in Get_zero.

**Integration_0_inf**: Do integration from 0 to infinity.

**Fourier_coefficients n**: Get the first n coefficients of Fourier series.

**FFT k il h**: where FFT means Fast Fourier Transform. n=pow(2,k) is the number of sampling points. il indicates whether to calculate amplitudes and angles, and h is the location of starting sample point (0.1*(i+h)). 

**random r u g n**: where random means we are going to create random values. r is the random seed. u is the average value. g2 is the variance of normal distribution, when g<=0, we ignore normal distribution restriction and just give random numbers in [0,1]. n is the number of data. 

**value a b h**: value means we just calculate the value of formula from a to b with step h.

## Class Calculator:
The structure of Calculator is listed below.
```cpp
class Calculator
{
public:
    Calculator(const string &expr,const string &deri); // initializing Calculator with strings.
    Calculator(const Calculator &value); // copy constructor. 
    
    void print();   // print expr
    void update(); // update the length of expr.
    bool num_detect(string str); // detect the type of input.
    void getVal(double &res); // res is the value of expr at specific point.
    double cal_num(const string &, double); // get an arithmetic value of expr.
    
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
    string::size_type pos, length;   // position for reading data and its length. 
    
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
```

## Main functions:
### double cal_num (const string &cal_expr, double x):
This function is used to get the numerical value of formula “cal_expr” at position x. This function calls getVal() to achieve its goal while getVal calls other functions in the private part of Caluculator to compute the value. The value would be returned as a double one.

**Algorithms:**  
The underlying mechanism of reading expressions from text is based on an operand stack and an operator stack. It calls readToken() to read and classify each word from input.txt, calls comparePrece() to compare the priority of different operators, then substitute independent variable with specific value to do the calculation.  

**Examples:**  
![GitHub Logo](/pictures/04.png) 
![GitHub Logo](/pictures/05.png)

### double integra (string& expr_i, double a, double b, double eps):
This function is used to calculate the integration of formula expr_i from a to b with precison eps. The result of integration will be returned as a double one.

**Algorithms:**  
The integration is based on the trapezoidal rule and its principle is clearly listed on Wiki Link:[trapezoidal rule](https://en.wikipedia.org/wiki/Trapezoidal_rule).

**Examples:**  
![GitHub Logo](/pictures/02.png) 
![GitHub Logo](/pictures/03.png)

### double integra_0i(string &expr_0i) 
This function is used to calculate the integration of formula expr_0i from 0 to infinity, and then return its result as a double value.

**Algorithms:**  
This integration is based on Gauss–Laguerre quadrature, which is an extension of Gaussian quadrature method for approximating the value of integrals of the following kind:  
![GitHub Logo](/pictures/06.png)  
This method is clearly listed on Wiki Link:[Gauss–Laguerre quadrature](https://en.wikipedia.org/wiki/Gauss%E2%80%93Laguerre_quadrature). We set the order of  Laguerre polynomial Ln(x) to be 5.

**Examples:**  
![GitHub Logo](/pictures/07.png) 
![GitHub Logo](/pictures/08.png)

### vector<double> zero(string &expr0,double a, double b, double h, double eps):
This function is used to find all the zero values of formula expr0 from a to b with step h. if f(x)<eps, then this x would be detected as a zero point. After finding all the zero values, all points will be stored in a vector as the return value.
                                                                                                         
**Algorithms:**                                                                                                    
This function is based on Bisection method, which is clearly listed on Wiki Link:[Bisection method](https://en.wikipedia.org/wiki/Bisection_method).
                                                                                                         
**Examples:**  
![GitHub Logo](/pictures/011.png) 
![GitHub Logo](/pictures/010.png)

### void limit1(double a, double b, double h, double eps):
This function directly read expression “orig and deri” from class and then find all the limit values of formula “orig” from a to b with step h. if f’(x)<eps, then this x would be detected as a zero point of f’(x), which is deri (the derivative of orig) in class. After finding all the zero values, all points will be printed out to the screen.
  
**Algorithms:**
For f(x), if f’(x)=0, then x is a point in which f(x) will have a limit maximum or minimum value. If there is no f’(x)=0 in the range from a to b, then local maximum and minimum value would be located on boundaries a and b. 

**Examples:**  
![GitHub Logo](/pictures/012.png) 
![GitHub Logo](/pictures/013.png)
![GitHub Logo](/pictures/014.png)

### void fourier(vector<double>::size_type n) 
This function is organized to get first n coefficients of Fourier series and then stored all results in a text file named “fourier_coefficients.txt”, where a and b are an and bn illustrated below.  
  
**Algorithms:**  
The introduction of Fourier series is clearly listed on Wiki Link:[Fourier series](https://en.wikipedia.org/wiki/Fourier_series).

**Examples:**    
![GitHub Logo](/pictures/015.png)
![GitHub Logo](/pictures/016.png)

### void fft(vector<double>::size_type k,vector<double>::size_type il,double h)
fft is used to do the Fast Fourier Transform. n=2k is the number of sample points, il = 0 means not to calculate magnitude and angle while = 1 means to calculate them, and h is the sample starting points. The result will be stored in a text file named “FFT_result.txt”.  
  
**Algorithms:**  
Please see Wiki Link:[FFT](https://en.wikipedia.org/wiki/Fast_Fourier_transform) to check how to do FFT.
  
**Examples:**    
![GitHub Logo](/pictures/017.png)
![GitHub Logo](/pictures/018.png)

### void random (double r,double u, double g, vector<double>::size_type n)
random is used to generate random values with designated average u and the variance g*g of normal distribution. r is the random seed. g^2 is the variance of normal distribution, when g<=0, we ignore normal distribution restriction and just create random numbers in the range [0,1]. n is the number of random data. The result will be stored in a text file named “random.txt”.

**Algorithms:**  
For generating random values between [0,1], m =2^16 and pi is the ith random value:  
![GitHub Logo](/pictures/022.png)  
For generating a bunch of random numbers that has average u and variance σ=g^2 under normal distribution, let n=12, then we get the lower formula. rndi are random numbers in [0,1]: 
![GitHub Logo](/pictures/023.png)
![GitHub Logo](/pictures/024.png)
  
**Examples:**    
![GitHub Logo](/pictures/020.png)
![GitHub Logo](/pictures/021.png)

### void value(double a, double b,double h)
This function is used to get the numerical value of orig in class Calculator.

**Examples:**    
![GitHub Logo](/pictures/01.png)
![GitHub Logo](/pictures/025.png)

### void cal_expression()
This function is not a class function but a one encapsulated in calculator.h. It is used to read input.txt and do corresponding calculation. The purpose of providing such a function is to simplify the users so that they can only write one code in main code and then get the result they need.

## Cautions:

1. This code does not distinguish minus ‘-’ and negative ‘-’. Although the code can successfully detect negative numbers, when there is input like “–exp(X)”, it will have an error.

2. The code take ‘X’ instead of ‘x’ as independent variable.

3. There is a limit on the value of f(x), when f(x)<1e-15 or f(x)>1e15, the code will throw an error. The range of limit could be modified by users themselves.
  
4. For input.txt, never change the sentences before ‘:’. You are only allowed to modify expressions and select specific operation you need.
  
5. All functions should be the same as that in cmath.
  
6. For integra_0i(), it only works for expressions that have the form as showed in its algorithms.
  
7. When using zero(), you must take step h seriously since it will affect the results.




