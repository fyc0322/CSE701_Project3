# CSE701_Project3 - Mathematical Calculator

## Author name: 

Yucheng Feng

## For users:

If you just want to use this program as a calculator to do mathematical calculation, read "Input" section and "Cautions" section. Then you will know how to use this one. For detailed information, read other parts of this doc and check comments in the code.

## Introduction:
“Mathematical Calculator” is fulfilled with C++ language and designed for mathematical operation in production-level. It can read formula from “input.txt” and then do sin, cos, tan, acos, asin, atan, cosh, sinh, tanh, acosh, asinh, atanh, log, log2, log10, exp, exp2, sqrt, pow, cbrt, hypot, fmax, fmin, abs, which are quite popular functions in C++ library “cmath”. It can also calculate integration, find limit-value, get zero-value, get Fourier coefficients, do FFT and create random points.

Codes for “Mathematical Calculator” are written in two files – calculator.h and Math_co.cpp. Most jobs are done in calculator.h, in which corresponding functions are encapsulated in Class “Calculator”. Math_co.cpp is just used to call a function to read formulas, print out or save calculation results.

“Mathematical Calculator” works for arbitrary formulas, but it has its own limitations which will be described in “Cautions” part of this document. Users of this code just need to give their formula in input.txt file, then run Math_co.cpp to get the results they are looking forward to.

## Input:
'#' means this line will be ignored.

 "/***** Input expressions *****/" means you should give your formulas after this.
 
 "/***** Operations *****/" means followings are the operations you want to do.

All input information should be listed in “input.txt”. This file has two parts. The first part is designed to read the formula and the second is used to do specific calculation. To read formulas from input file, users should give their own mathematical expressions. The format of basic functions should be in the same format as those in the C++ cmath library. The second is to select specific operations. Nowadays, this code offers eight possible operations; users should choose their own operations based on their needs. The format of each operation is listed below:

**Integration a b eps**: where Integration indicates that we are going to do integration from a to b with precision eps.

**Get_zero a b h eps**: where Get_zero indicates that we are going to find all the zero values of input formula in the range from a to b with step h and precision eps. Get_zero is based on Bisection method so it needs a step.

**Limit a b  h  eps**: where Limit1 indicates that we are going to find the limit value of formula in the given range and other input parameters’ meaning are just like these in Get_zero.

**Integration_0_inf**: Do integration from 0 to infinity.

**Fourier_coefficients n**: Get the first n coefficients of Fourier series.

**FFT k il h**: where FFT means Fast Fourier Transform. n=pow(2,k) is the number of sampling points. il indicates whether to calculate amplitudes and angles, and h is the location of starting sample point (0.1*(i+h)). 

**Random u g n**: where random means we are going to create random values. u is the average value. g2 is the variance of normal distribution, when g<=0, we ignore normal distribution restriction and just give random numbers in [0,1]. n is the number of data. 

**Value a b h**: value means we just calculate the value of formula from a to b with step h.

## Class Calculator:
The structure of Calculator is listed below.
```cpp
/**
 * @brief class with strings to store formulas and functions to do mathematical calculations. 
 */
class Calculator
{
public:
    Calculator(const string &expr, const string &deri);                            // Initializing Calculator with a string.
    Calculator(const Calculator &value);                                           // Copy constructor to initialize an object of Calculator with another existing Calculator object "value".
    void print();                                                                  // Print mathematical expressions out.
    void update();                                                                 // Update the length of expr.
    bool num_detect(string str);                                                   // Detect the type of input.
    void getVal(double &res);                                                      // Calculate the value of res.
    double cal_num(const string &cal_expr, double x);                              // Get an function value of cal_expr at point x.
    double integra(string &expr_i, double a, double b, double eps);                // Integration operation.
    double integra_0i(string &expr_0i);                                            // Do integration from 0 to infinity.
    vector<double> zero(string &expr0, double a, double b, double h, double eps);  // Get zero value.
    void limit(double a, double b, double h, double eps);                          // Get local limit value with f(x) and f'(x).
    void fourier(vector<double>::size_type n);                                     // Get fourier coefficients.
    void fft(vector<double>::size_type k, vector<double>::size_type il, double h); // Fast Fourier Transform.
    void random(double u, double g, vector<double>::size_type n);                  // Create n random values with average u and variance g.
    void value(double a, double b, double h);                                      // Calculate function values.

private:
/**
 * @param expr Store the mathematical expression for calculation.
 * @param orig Store original mathematical expression.
 * @param deri Store the derivative expression of orig.
 * @param token The word read at each time.
 * @param tkType Word type. 
 * @param pos Position for read. 
 * @param length The length of formula.
 * @param ptrList List of supported operators.
 * @param ptrArgCnt The number of arguments required for the operator.
 * @param preceMap Operation priority table.
 * 
 */
    string expr;
    string orig;
    string deri;
    string token;
    TokenType tkType;
    string::size_type pos, length;
    static string ptrList[];
    static int32_t ptrArgCnt[];
    static int32_t preceMap[][ptr_num];

    /**
     * @brief Read next word.
     */
    void readToken();

    /**
     * @brief Compare the priority of the two operators (ptr1 and ptr2), and then return an int32_t value.
     * 
     * @param ptr1 Operator 1.
     * @param ptr2 Operator 2.
     * @return int32_t Return the result of priority.
     */
    int32_t comparePrece(const string &ptr1, const string &ptr2);

    /**
     * @brief Single-step calculation of an operator and return a double value.
     * 
     * @param ptr The operator.
     * @param arg Arguments required by the operator.
     * @return double Value after calculation. 
     */
    double calculate(const string &ptr, double arg[]);

    /**
     * @brief Single-step calculation of a function and return a double value.
     * 
     * @param fun The function.
     * @param arg Arguments required by the function.
     * @return double Value after calculation. 
     */
    double callFun(const string &fun, double arg[]);

    /**
     * @brief Get operator (ptr) sequence number, and then return an int32_t value.
     * 
     * @param ptr The operator.
     * @return int32_t The operator's sequence number.
     */
    int32_t getPtrIndex(const string &ptr);

    /**
     * @brief Detect whether the function is supported or not.
     * 
     * @param fun The function.
     * @return An object of enum.
     */
    funList getFunIndex(const string &fun);

    /**
     * @brief Check whether the number of arguments in opnd matches n. If match, get n parameters from operand stack (opnd) and store them in arg.
     * 
     * @param opnd The operator/function.
     * @param arg The result.
     * @param n Number of arguments needed to do a calculation.
     * @return true The number of arguments in opnd matches n.
     * @return false The number of arguments in opnd fails to match n.
     */
    bool getArg(stack<double> &opnd, double arg[], int32_t n);
};
```

## Main functions:
### double cal_num (const string &cal_expr, double x):
This function is used to get the numerical value of formula “cal_expr” at position x. This function calls getVal() to achieve its goal while getVal calls other functions in the private part of Calculator to compute the value. The value would be returned as a double one.

**Algorithms:**  
The underlying mechanism of reading expressions from text is based on an operand stack and an operator stack. It calls readToken() to read and classify each word from input.txt, calls comparePrece() to compare the priority of different operators, then substitute independent variable with specific value to do the calculation.  

**Examples:**  
```
f(X)=-exp(-1+2)   
f(X)=-2.71828
```
```
f(X)=-exp(-1+2)+fmax(1,3)+log(10)*sin(3)   
f(X)=0.606659
```

### double integra (string& expr_i, double a, double b, double eps):
This function is used to calculate the integration of formula expr_i from a to b with precision eps. The result of integration will be returned as a double one.

**Algorithms:**  
The integration is based on the trapezoidal rule and its principle is clearly listed on Wiki Link:[trapezoidal rule](https://en.wikipedia.org/wiki/Trapezoidal_rule).

**Examples:**  
```
f(X)=exp(-X)
***** Integration Operation *****
integration from 0 to 1 of f(X) is: 0.632121
```
```
f(X)=X*X
***** Integration Operation *****
integration from 0 to 1 of f(X) is: 0.333333
```

### double integra_0i(string &expr_0i) 
This function is used to calculate the integration of formula expr_0i from 0 to infinity, and then return its result as a double value.

**Algorithms:**  
This integration is based on Gauss–Laguerre quadrature, which is an extension of Gaussian quadrature method for approximating the value of integrals of the following kind:  
![GitHub Logo](/pictures/06.png)  
This method is clearly listed on Wiki Link:[Gauss–Laguerre quadrature](https://en.wikipedia.org/wiki/Gauss%E2%80%93Laguerre_quadrature). We set the order of  Laguerre polynomial Ln(x) to be 5.

**Examples:**  
```
f(X)=X*exp(-X)
***** Integrate from 0 to infinity *****
integration from 0 to  infinity of f(X) is: 1
```
```
f(X)=sin(X)*exp(-X)
***** Integrate from 0 to infinity *****
integration from 0 to  infinity of f(X) is: 0.498905
```

### vector<double> zero(string &expr0,double a, double b, double h, double eps):
This function is used to find all the zero values of formula expr0 from a to b with step h. if f(x)<eps, then this x would be detected as a zero point. After finding all the zero values, all points will be stored in a vector as the return value.
                                                                                                         
**Algorithms:**                                                                                                    
This function is based on Bisection method, which is clearly listed on Wiki Link:[Bisection method](https://en.wikipedia.org/wiki/Bisection_method).
                                                                                                         
**Examples:**  
```
f(X)=2-exp(X*X)
***** Get zero Operation *****
x-points for 0 from -1 to 1 are: -0.832555 0.832555 
```
```
f(X)=2-exp(-X*X)
***** Get zero Operation *****
function at interval -1 to 1 has no zero value.
```

### void limit(double a, double b, double h, double eps):
This function directly read expression “orig and deri” from class and then find all the limit values of formula “orig” from a to b with step h. if f’(x)<eps, then this x would be detected as a zero point of f’(x), which is deri (the derivative of orig) in class. After finding all the zero values, all points will be printed out to the screen.
  
**Algorithms:**
For f(x), if f’(x)=0, then x is a point in which f(x) will have a limit maximum or minimum value. If there is no f’(x)=0 in the range from a to b, then local maximum and minimum value would be located on boundaries a and b. 

**Examples:**  
```
f(X)=2+exp(-X*X)
f'(X)=-2*X*exp(-X*X)
***** Limit Operation *****
x-points for 0 from -1 to 1 are: 0 
0 is the local maxima point. The maxima is: 3
````
```
f(X)=X*X
f'(X)=2*X
***** Limit Operation *****
x-points for 0 from -1 to 1 are: 0 
0 is the local minimum point. The minimum is: 0
````
```
f(X)=exp(X)
f'(X)=exp(X)
***** Limit Operation *****
function at interval -1 to 1 has no zero value.
no limit value between -1 to 1. largest and lowest value are both on the boundaries
highest value is: 2.71828 lowest value is: 0.367879
````
    
### void fourier(vector<double>::size_type n) 
This function is organized to get first n coefficients of Fourier series and then stored all results in a text file named “fourier_coefficients.txt”, where a and b are an and bn illustrated below.  
  
**Algorithms:**  
The introduction of Fourier series is clearly listed on Wiki Link:[Fourier series](https://en.wikipedia.org/wiki/Fourier_series).

**Examples:**    
```
f(X)=X*X
***** Fourier_coefficients Operation *****
Successfully calculate Fourier_coefficients and store results in fourier_coefficients.txt.
For some results in fourier_coefficients.txt:
 i         a              b
 0          26.3183       0
 1          3.60656       -12.6867
 2         0.606573       -6.33721
 3        0.0509917       -4.21799
 4        -0.143439       -3.15631
 5        -0.233435       -2.51767
 6        -0.282321       -2.09052
 7        -0.311775       -1.78423
 8        -0.330917       -1.55347
 9        -0.343983       -1.37302
10        -0.353376       -1.22808
````
````
f(X)=exp(X)
***** Fourier_coefficients Operation *****
Successfully calculate Fourier_coefficients and store results in fourier_coefficients.txt.
For some results in fourier_coefficients.txt:
 i         a              b
 0          170.107       0
 1          82.3527       -87.698
 2             29.7       -70.0906
 3           12.149       -52.4831
 4          4.92224       -41.0698
 5          1.34866       -33.4685
 6         -0.65801       -28.1207
 7         -1.89092       -24.1741
 8         -2.70059       -21.1468
 9         -3.25953       -18.7521
10         -3.66281       -16.8103
````

### void fft(vector<double>::size_type k,vector<double>::size_type il,double h)
fft is used to do the Fast Fourier Transform. n=2k is the number of sample points, il = 0 means not to calculate magnitude and angle while = 1 means to calculate them, and h is the sample starting points. The result will be stored in a text file named “FFT_result.txt”.  
  
**Algorithms:**  
Please see Wiki Link:[FFT](https://en.wikipedia.org/wiki/Fast_Fourier_transform) to check how to do FFT.
  
**Examples:**    
````
f(X)=exp(X)
***** FFT Operation *****
Successfully calculate FFT and store results in FFT_result.txt.
For some results in FFT_result.txt:
 f(X)       "Real part after FFT"   Imag part after FFT 
   1.05127        6005.95                     0
   1.16183        2905.87               3152.55
   1.28403        990.395               2544.07
   1.41907        342.366               1907.45
   1.56831        74.1327               1489.82
   1.73325       -58.8019               1209.58
   1.91554       -133.526               1011.11
     2.117       -179.463               863.634
   2.33965       -209.631               749.662
   2.58571       -230.466               658.741
````
````
f(X)=exp(-X)
***** FFT Operation *****
Successfully calculate FFT and store results in FFT_result.txt.
For some results in FFT_result.txt:
 f(X)       "Real part after FFT"   Imag part after FFT 
  0.951229        9.97923                     0
  0.860708        5.31844              -4.73967
  0.778801        2.43865              -3.82485
  0.704688        1.46438              -2.86773
  0.637628         1.0611              -2.23986
   0.57695       0.861244              -1.81854
  0.522046         0.7489              -1.52015
  0.472367       0.679837              -1.29842
  0.427415       0.634482              -1.12707
````

### void random (double u, double g, vector<double>::size_type n)
random is used to generate random values with designated average u and the variance g*g of normal distribution. g^2 is the variance of normal distribution, when g<=0, we ignore normal distribution restriction and just create random numbers in the range [0,1]. n is the number of random data. The result will be stored in a text file named “random.txt”.

**Algorithms:**  
For generating random values between [0,1], m =2^16 and pi is the ith random value:  
<img src="/pictures/022.png" width="400" height="60">  
For generating a bunch of random numbers that has average u and variance σ=g^2 under normal distribution, let n=12, then we get the lower formula. rndi are random numbers in [0,1]:   
<img src="/pictures/023.png" width="300" height="100">   
<img src="/pictures/024.png" width="300" height="60">  
  
**Examples:**    
````
***** random Operation *****
random points in [0,1] are: 
0.592117,  0.828156,  0.414597,  0.378052,  0.351578,  0.000457764,  0.151108,  0.435608,  0.514359,  0.189362,  0.970627,  0.908203,  0.752335,  0.754242,
0.670029,  0.780212,  0.987381,  0.304474,  0.296188,  0.286011,  0.391373,  0.699432,  0.145981,  0.909973,  0.386185,  0.0484924,  0.766281,  0.386475, 
0.643692,  0.711029,  0.953964,  0.69989,  0.0857697,  0.296417,  0.755905,  0.0845947,  0.884293,  0.664032,  0.468979,  0.0249634,  0.461136,  0.923248, 
0.64006,  0.255371,  0.488174,  0.433441,  0.0660248,  0.760193,  0.887283,  0.803986,  End
````
````
***** random Operation *****
Normal distribution random points with average= 1 variance= 1.5:
0.418411,  2.36848,  1.51192,  1.31746,  2.25386,  3.28987,  0.394241,  -1.46428,  1.18306,  1.30501,  -0.629684,  1.84773,  1.70601,  0.913895,  -0.0598602, 
0.753494,  -0.677292,  -0.883469,  -0.896286,  -0.246994,  0.0331573,  -1.08708,  1.36104,  3.34627,  0.837357,  0.303055,  -0.787888,  -0.466721,  -1.26469,
-1.21306,  0.156937,  3.31404,  1.22701,  1.86458,  -0.304489,  1.18855,  3.81245,  0.535965,  0.827835,  0.656815,  0.491653,  0.801102,  2.05391,  3.21883,
-1.2354,  1.15999,  0.373734,  -3.12541,  0.131302,  1.61263,  End
````

### void value(double a, double b,double h)
This function is used to get the numerical value of orig in class Calculator.

**Examples:**    
````
f(X)=-exp(X)+1
***** value Operation *****
f(0)= 0
f(1)= -1.71828
f(2)= -6.38906
f(3)= -19.0855
````
````
f(X)=-exp(X)+cos(X)-fmin(2,4)
***** value Operation *****
f(0)= -2
f(1)= -4.17798
f(2)= -9.8052
f(3)= -23.0755
````

### void cal_expression()
This function is not a class function but a one encapsulated in calculator.h. It is used to read input.txt and do corresponding calculation. The purpose of providing such a function is to simplify the users so that they can only write one code in main code and then get the result they need.

## Cautions:

1. The code take ‘X’ instead of ‘x’ as independent variable.

2. There is a limit on the value of f(x), when f(x)<1e-15 or f(x)>1e15, the code will throw an invalid input error. The range of the limit could be modified by users themselves.
  
3. All functions should be in the same format as these listed in cmath.
  
4. For integra_0i(), it only works for expressions that have the form as showed in its algorithms.
  
5. When using functions concerning step h, be careful that currently the smallest step is set to be 0.001. Users could change this restriction.



