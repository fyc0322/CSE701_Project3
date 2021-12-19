/**
 * @file calculator.h
 * @author Yucheng Feng (fengy126@mcmaster.ca)
 * @brief Mathematical Calculator
 * @version 0.1
 * @date 2021-12-06
 * 
 * @copyright Copyright (c) 2021
 */
#include <string>
#include <stack>
#include <cmath>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <ctime>
#include <map>
using namespace std;

/**
 * @param step_limit The minimal step length should be 10^-3.
 * @param ptr_num Number of supported operators.
 * @param fun_num Number of supported functions.
 */
const int32_t step_limit=3;
const string::size_type ptr_num = 9;
const string::size_type fun_num = 24;

/**
 * @brief TokenType gathers types of words in the expression.
 * @param TKT_NUMBER Number.
 * @param TKT_OPERATOR Operator.
 * @param  TKT_FUNCTION Function.
 * @param TKT_ENDSIGN Terminator.
 * @param TKT_UNKNOW Unknown symbol.
 */
typedef enum
{
    TKT_NUMBER,
    TKT_OPERATOR,
    TKT_FUNCTION,
    TKT_ENDSIGN,
    TKT_UNKNOW
} TokenType;

/**
 * @brief funList collects all supported functions.
 */
enum funList
{
    sina = 1,
    cosa,
    tana,
    acosa,
    asina,
    atana,
    cosha,
    sinha,
    tanha,
    acosha,
    asinha,
    atanha,
    loga,
    log2a,
    log10a,
    expa,
    exp2a,
    sqrta,
    powa,
    cbrta,
    hypota,
    fmaxa,
    fmina,
    absa,
};

/**
 * @brief fun_str maps a string to its corresponding function.
 */
map<string, funList> fun_str{
    {"sin(", sina},
    {"cos(", cosa},
    {"tan(", tana},
    {"acos(", acosa},
    {"asin(", asina},
    {"atan(", atana},
    {"cosh(", cosha},
    {"sinh(", sinha},
    {"tanh(", tanha},
    {"acosh(", acosha},
    {"asinh(", asinha},
    {"atanh(", atanha},
    {"log(", loga},
    {"log2(", log2a},
    {"log10(", log10a},
    {"exp(", expa},
    {"exp2(", exp2a},
    {"sqrt(", sqrta},
    {"pow(", powa},
    {"cbrt(", cbrta},
    {"hypot(", hypota},
    {"fmax(", fmaxa},
    {"fmin(", fmina},
    {"abs(", absa},
};

/**
 * @brief funArgCnt maps a function with how many parameters it needs.
 */
map<string, int32_t> funArgCnt{
    {"sin(", 1},
    {"cos(", 1},
    {"tan(", 1},
    {"acos(", 1},
    {"asin(", 1},
    {"atan(", 1},
    {"cosh(", 1},
    {"sinh(", 1},
    {"tanh(", 1},
    {"acosh(", 1},
    {"asinh(", 1},
    {"atanh(", 1},
    {"log(", 1},
    {"log2(", 1},
    {"log10(", 1},
    {"exp(", 1},
    {"exp2(", 1},
    {"sqrt(", 1},
    {"pow(", 2},
    {"cbrt(", 1},
    {"hypot(", 2},
    {"fmax(", 2},
    {"fmin(", 2},
    {"abs(", 1},
};

/**
 * @brief Class with strings to store formulas and functions to do mathematical calculations. 
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

/**
 * @brief Initialize static data. "f(" means function and "#" is the terminator.
 */
string Calculator::ptrList[] = {"f(", ",", "+", "-", "*", "/", "(", ")", "#"};

int32_t Calculator::ptrArgCnt[] = {0, 1, 2, 2, 2, 2, 0, 0, 0};

int32_t Calculator::preceMap[][ptr_num] =
    {
        //{f(, , +, -, *, /, (, ), #}
        {-1, -1, -1, -1, -1, -1, -1, 1, 1}, // f(
        {1, 1, 1, 1, 1, 1, 1, 1, 1},        // ,
        {-1, 1, 1, 1, -1, -1, -1, 1, 1},    // +
        {-1, 1, 1, 1, -1, -1, -1, 1, 1},    // -
        {-1, 1, 1, 1, 1, 1, -1, 1, 1},      // *
        {-1, 1, 1, 1, 1, 1, -1, 1, 1},      // /
        {-1, 1, -1, -1, -1, -1, -1, 0, 1},  // (
        {1, 1, 1, 1, 1, 1, 1, 1, 1},        // )
        {-1, 1, -1, -1, -1, -1, -1, -1, 0}, // #
};

/**
 * @brief Initializing Calculator with strings.
 * 
 * @param expr Store the mathematical expression for calculation.
 * @param orig Store original mathematical expression.
 * @param deri Store the derivative expression of orig.
 * @param token  The word read at each time.
 * @param pos Position for reading data
 * @param length The length of formula
 */
Calculator::Calculator(const string &expr, const string &deri)
{
    this->expr = expr + "#";
    this->orig = expr + "#";
    this->deri = deri + "#";
    this->token = "";
    this->pos = 0;
    this->length = 1 + expr.length();
}

/**
 * @brief Construct a new Calculator:: Calculator object with an existing Calculator:: Calculator object "value".
 * 
 * @param value An existing Calculator:: Calculator object.
 */
Calculator::Calculator(const Calculator &value)
{
    this->expr = value.expr;
    this->orig = value.orig;
    this->deri = value.deri;
    this->token = value.token;
    this->pos = 0;
    this->length = value.length;
}

/**
 * @brief Print out orig and deri.
 */
void Calculator::print()
{
    cout << orig << endl;
    cout << deri << endl;
}

/**
 * @brief Update the length of expr.
 */
void Calculator::update()
{
    length = expr.length();
}

/**
 * @brief Detect whether fstr is a number or a letter.
 * 
 * @param fstr Input string.
 * @return true fstr is a number.
 * @return false fstr is not a number.
 */
bool Calculator::num_detect(string fstr)
{
    stringstream str_(fstr);
    double d;
    string c;
    if (!(str_ >> d))
        return false;
    if (str_ >> c)
        return false;
    return true;
}

/**
 * @brief Get operator (ptr) sequence number.
 * 
 * @param ptr The operator.
 * @return int32_t Sequence number of the operator (ptr).
 */
int32_t Calculator::getPtrIndex(const string &ptr)
{
    for (string::size_type i = 0; i < ptr_num; ++i)
    {
        if (ptrList[i] == ptr)
            return static_cast<int32_t>(i);
    }
    return -1;
}

/**
 * @brief Get function (fun) sequence number.
 * 
 * @param fun The function.
 * @return int32_t Sequence number of the operator (ptr)
 */
funList Calculator::getFunIndex(const string &fun)
{
    funList t = fun_str[fun];
    return t;
}

/**
 * @brief Compare the priority of the two operators (ptr1 and ptr2).
 * 
 * @param ptr1 Operator 1.
 * @param ptr2 Operator 2.
 * @return int32_t ptr1's priority compared with ptr2.
 */
int32_t Calculator::comparePrece(const string &ptr1, const string &ptr2)
{
    int32_t m = getPtrIndex(ptr1);
    int32_t n = getPtrIndex(ptr2);
    // if m or n == -1, then m or n is a function.
    if (m == -1)
        m = 0;
    if (n == -1)
        n = 0;
    return preceMap[m][n];
}

/**
 * @brief Single-step calculation of an operator.
 * 
 * @param ptr An operator.
 * @param arg Arguments for the operator.
 * @return double Result of calculation.
 */
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

/**
 * @brief Single-step calculation of a function and set value limit.
 * 
 * @param fun The function.
 * @param arg Arguments for the function.
 * @return double Result of such a function.
 */
double Calculator::callFun(const string &fun, double arg[])
{
    switch (getFunIndex(fun))
    {
    case sina:
        return sin(arg[0]);

    case cosa:
        return cos(arg[0]);

    case tana:
        return tan(arg[0]);

    case acosa:
        if (arg[0] > 1 || arg[0] < -1)
        {
            cout << "Error: invalid input for acos().\n"
                 << endl;
            exit(3);
        }
        return acos(arg[0]);

    case asina:
        if (arg[0] > 1 || arg[0] < -1)
        {
            cout << "Error: invalid input for asin().\n"
                 << endl;
            exit(3);
        }
        return asin(arg[0]);

    case atana:
        return atan(arg[0]);

    case cosha:
        if (cosh(arg[0]) > 1e15)
        {
            cout << "Error: invalid input for cosh().\n"
                 << endl;
            exit(3);
        }
        return cosh(arg[0]);

    case sinha:
        if (sinh(arg[0]) > 1e15)
        {
            cout << "Error: invalid input for sinh().\n"
                 << endl;
            exit(3);
        }
        return sinh(arg[0]);

    case tanha:
        if (tanh(arg[0]) < 1e-15)
        {
            cout << "Error: invalid input for tanh().\n"
                 << endl;
            exit(3);
        }
        return tanh(arg[0]);

    case acosha:
        if (arg[0] < 1)
        {
            cout << "Error: invalid input for acosh().\n"
                 << endl;
            exit(3);
        }
        return acosh(arg[0]);

    case asinha:
        return asinh(arg[0]);

    case atanha:
        if (arg[0] >= 1 || arg[0] <= -1)
        {
            cout << "Error: invalid input for atanh().\n"
                 << endl;
            exit(3);
        }
        return atanh(arg[0]);

    case loga:
        if (arg[0] <= 0)
        {
            cout << "Error: invalid input for log().\n"
                 << endl;
            exit(3);
        }
        return log(arg[0]);

    case log2a:
        if (arg[0] <= 0)
        {
            cout << "Error: invalid input for log2().\n"
                 << endl;
            exit(3);
        }
        return log2(arg[0]);

    case log10a:
        if (arg[0] <= 0)
        {
            cout << "Error: invalid input for log10().\n"
                 << endl;
            exit(3);
        }
        return log10(arg[0]);

    case expa:
        if (exp(arg[0]) > 1e15 || exp(arg[0]) < 1e-15)
        {
            cout << "Error: invalid input for exp().\n"
                 << endl;
            exit(3);
        }
        return exp(arg[0]);

    case exp2a:
        if (exp2(arg[0]) > 1e15 || exp2(arg[0]) < 1e-15)
        {
            cout << "Error: invalid input for exp2().\n"
                 << endl;
            exit(3);
        }
        return exp2(arg[0]);

    case sqrta:
        if (arg[0] < 0)
        {
            cout << "Error: invalid input for sqrt().\n"
                 << endl;
            exit(3);
        }
        return sqrt(arg[0]);

    case powa:
        return pow(arg[0], arg[1]);

    case cbrta:
        return cbrt(arg[0]);

    case hypota:
        if (arg[0] <= 0 || arg[1] <= 0)
        {
            cout << "Error: invalid input for hypot().\n"
                 << endl;
            exit(3);
        }
        return hypot(arg[0], arg[1]);

    case fmaxa:
        return fmax(arg[0], arg[1]);

    case fmina:
        return fmin(arg[0], arg[1]);

    case absa:
        return abs(arg[0]);
    }
    return 0;
}

/**
 * @brief Read next word from input.txt
 */
void Calculator::readToken()
{
    // Check whether the reading is completed.
    if (pos >= length)
    {
        tkType = TKT_ENDSIGN;
        return;
    }
    string::size_type pos_t = pos;
    char ch = expr[pos_t++];
    // Judge the type of ch.
    if (-1 != getPtrIndex(std::string(1, ch)))
    {
        if (ch != '-')
        {
            tkType = TKT_OPERATOR;
        }
        else
        {
            if (pos_t > 1 && (expr[pos_t - 2] == ')' || isdigit(expr[pos_t - 2])))
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

    else if (isdigit(ch))
    {
        while (pos_t < length && isdigit(ch = expr[pos_t]))
        {
            ++pos_t;
        }

        if (ch == '.')
        {
            ++pos_t;
            while (pos_t < length && isdigit(expr[pos_t]))
                ++pos_t;
        }
        tkType = TKT_NUMBER;
    }

    else if (isalpha(ch))
    {
        while (pos_t < length && (isalnum(ch) || ch == '_'))
        {
            ch = expr[++pos_t];
        }

        if (ch == '(')
        {
            ++pos_t;
            if (getFunIndex(expr.substr(pos, pos_t - pos)))
                tkType = TKT_FUNCTION;
            else
                tkType = TKT_UNKNOW;
        }
        else
            tkType = TKT_UNKNOW;
    }
    else
        tkType = TKT_UNKNOW;
    token = expr.substr(pos, pos_t - pos);
    pos = pos_t;
}

/**
 * @brief Check whether the number of arguments in opnd matches n. If match, get n parameters from operand stack (opnd) and store them in arg.
 * 
 * @param opnd The operand stack.
 * @param arg Store required number of operands.
 * @param n Required number of operands.
 * @return true The number of arguments in opnd matches n.
 * @return false The number of arguments in opnd does not match n.
 */
bool Calculator::getArg(stack<double> &opnd, double arg[], int32_t n)
{
    if (opnd.size() < static_cast<unsigned long>(n))
    {
        return false;
    }

    for (int32_t i = n - 1; i >= 0; --i)
    {
        arg[i] = opnd.top();
        opnd.pop();
    }
    return true;
}

/**
 * @brief Solve expr and return its value in res.
 * @param res Return value of function value.
 * @param optr Stack for operators.
 * @param opnd Stack for operands.
 * @param comRes Store priority of two operators.
 * @param idx Store priority of two operators.
 * @param argCnt Store the required number of operands to operate a function or an operator.
 */
void Calculator::getVal(double &res)
{
    stack<string> optr;
    stack<double> opnd;
    int32_t comRes;
    int32_t idx;
    int32_t argCnt;
    optr.push("#");
    pos = 0;
    readToken();
    // Calculate the value of input expression.
    while (tkType != TKT_ENDSIGN || !optr.empty())
    {

        if (tkType == TKT_UNKNOW)
        {
            cout << "Error: The input mathematical expression contains unknown characters\n"
                 << endl;
            exit(1);
        }
        // If the word from readToken is a number, then push it in opnd.
        if (tkType == TKT_NUMBER)
        {
            opnd.push(atof(token.c_str()));
            readToken();
        }
        // Do the calculation based on operation priority.
        else
        {
            comRes = comparePrece(optr.top(), token);
            switch (comRes)
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
                            cout << "Error: The number of operands failed to match the requirement\n"
                                 << endl;
                            exit(2);
                        }
                        res = calculate(ptr, arg);
                        opnd.push(res);
                    }
                }
                else
                {
                    idx = getFunIndex(ptr);
                    argCnt = funArgCnt[ptr];
                    if (!getArg(opnd, arg, argCnt))
                    {
                        cout << "Error: The number of operands failed to match the requirement\n"
                             << endl;
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

/**
 * @brief Calculate the function value at given point x.
 * 
 * @param cal_expr Formula f(X) waiting for solving.
 * @param x Given point.
 * @return double The value of f(x).
 */
double Calculator::cal_num(const string &cal_expr, double x)
{
    double res;
    if (abs(x) < 1e-8)
    {
        x = 0;
    }
    stringstream ss;
    if (cal_expr.find('#') == string::npos)
    {
        expr = cal_expr + "#";
    }
    else
        expr = cal_expr;
    ss << setiosflags(ios::fixed) << setprecision(step_limit) << x;
    string str_n = ss.str();
    while (expr.find('X') < expr.length())
    {
        expr.replace(expr.find('X'), 1, str_n);
    }
    update();
    getVal(res);
    return res;
}

/**
 * @brief Do integration of expr_i from a to b with precision eps.
 * 
 * @param expr_i Formula waiting for integration.
 * @param a Lower boundary for integration.
 * @param b Upper boundary for integration.
 * @param eps Precision.
 * @return double The result of integration.
 */
double Calculator::integra(string &expr_i, double a, double b, double eps)
{
    // n: divide the integral interval into 2n equal parts.
    // fa and fb are the lower and upper boundary.
    // h is the width of the integral interval.
    // t1 is the return value of integration result.
    // other variables are for temporary usage.
    int32_t n, k;
    double fa, fb, h, t1, p, s, x, t = 10, l;
    if (a > b)
    {
        cout << "Integration Error: upper limit is lower than lower limit: " << endl;
        exit(2);
    }

    // Get function value at left and right boundary.
    fa = cal_num(expr_i, a);
    fb = cal_num(expr_i, b);
    // Do the integration.
    n = 1;
    h = b - a;
    t1 = h * (fa + fb) / 2.0;
    p = 1.0;
    while (p >= eps)
    {
        s = 0.0;
        for (k = 0; k < n; k++)
        {
            x = a + (k + 0.5) * h;
            l = cal_num(expr_i, x);
            s = s + l;
        }
        t = (t1 + h * s) / 2.0;
        p = fabs(t1 - t);
        t1 = t;
        n = n + n;
        h = h / 2.0;
    }
    return (t);
}

/**
 * @brief Do integration from 0 to infinity.
 * 
 * @param expr_0i Formula waiting for integration.
 * @return double The result of integration.
 */
double Calculator::integra_0i(string &expr_0i)
{
    // x are the roots of Laguerre polynomial Ln(x).
    // lc are the weights for integration.
    // g is the integration result .
    vector<double>::size_type i;
    double x, g;
    vector<double> t{0.26355990, 1.41340290, 3.59642600, 7.08580990, 12.64080000};
    vector<double> Lc{0.6790941054, 1.638487956, 2.769426772, 4.315944000, 7.104896230};
    g = 0.0;
    for (i = 0; i <= 4; i++)
    {
        x = t[i];
        g = g + Lc[i] * cal_num(expr_0i, x);
    }
    return (g);
}

/**
 * @brief Get the zero-value of expr0.
 * 
 * @param expr0 Formula waiting for solving zero-solution.
 * @param a Lower boundary for zero-solution.
 * @param b Upper boundary for zero-solution.
 * @param h Step length.
 * @param eps Precision. f(x)<eps indicates x is zero-point.
 * @return vector<double> Return all zero-points.
 */
vector<double> Calculator::zero(string &expr0, double a, double b, double h, double eps)
{
    // n is the number of zero-points.
    // x is a group of zero points.
    // y, y0 and y1 are the function values at end points.
    // z, z0 and z1 are  end points.
    int n, js;
    vector<double> x;
    double z, y, z1, y1, z0, y0;
    n = 0;
    z = a;
    y = cal_num(expr0, z);
    // Bisection method to get zero-points.
    while ((z <= b + h / 2.0))
    {
        if (fabs(y) < eps)
        {
            n = n + 1;
            x.push_back(z);
            z = z + h / 2.0;
            y = cal_num(expr0, z);
        }
        else
        {
            z1 = z + h;
            y1 = cal_num(expr0, z1);
            if (fabs(y1) < eps)
            {
                n = n + 1;
                x.push_back(z1);
                z = z1 + h / 2.0;
                y = cal_num(expr0, z);
            }
            else if (y * y1 > 0.0)
            {
                y = y1;
                z = z1;
            }
            else
            {
                js = 0;
                while (js == 0)
                {
                    if (fabs(z1 - z) < eps)
                    {
                        n = n + 1;
                        x.push_back((z1 + z) / 2.0);
                        z = z1 + h / 2.0;
                        y = cal_num(expr0, z);
                        js = 1;
                    }
                    else
                    {
                        z0 = (z1 + z) / 2.0;
                        y0 = cal_num(expr0, z0);
                        if (fabs(y0) < eps)
                        {
                            x.push_back(z0);
                            n = n + 1;
                            js = 1;
                            z = z0 + h / 2.0;
                            y = cal_num(expr0, z);
                        }
                        else if ((y * y0) < 0.0)
                        {
                            z1 = z0;
                            y1 = y0;
                        }
                        else
                        {
                            z = z0;
                            y = y0;
                        }
                    }
                }
            }
        }
    }
    for (vector<double>::size_type i = 0; i < x.size(); i++)
        if (abs(x[i]) < 1e-10)
        {
            x[i] = 0;
        }
    if (n == 0)
    {
        cout << "function at interval " << a << " to " << b << " has no zero value.\n"
             << endl;
    }
    if (n > 0)
    {
        cout << "x-points for 0 from " << a << " to " << b << " are: ";
        for (vector<double>::size_type i = 0; i < x.size(); i++)
        {
            cout << x[i] << " ";
        }
        cout << "\n"
             << endl;
    }
    return (x);
}

/**
 * @brief Find all the limit values of formula orig.
 * 
 * @param a Lower boundary for limit-solution.
 * @param b Upper boundary for limit-solution.
 * @param h Step length.
 * @param eps Precision.
 */
void Calculator::limit(double a, double b, double h, double eps)
{
    vector<double> dr_zero;
    vector<double>::size_type n;
    double lb, rb;
    dr_zero = zero(deri, a, b, h, eps);
    n = dr_zero.size();
    if (n == 0)
    {
        cout << "no limit value between " << a << " to " << b << ". largest and lowest value are both on the boundaries" << endl;
        lb = cal_num(orig, a);
        rb = cal_num(orig, b);
        (lb > rb) ? cout << "highest value is: " << lb << " lowest value is: " << rb << "\n"
                         << endl
                  : cout << "highest value is: " << rb << " lowest value is: " << lb << "\n"
                         << endl;
    }
    else
    {
        for (vector<double>::size_type i = 0; i < dr_zero.size(); i++)
        {
            lb = cal_num(orig, dr_zero[i]);
            rb = cal_num(orig, dr_zero[i] + 0.1);
            (lb > rb) ? cout << dr_zero[i] << " is the local maxima point. The maxima is: " << lb << endl : cout << dr_zero[i] << " is the local minimum point. The minimum is: " << lb << "\n"
                                                                                                                 << endl;
        }
    }
}

/**
 * @brief Get Fourier series coefficients and store them in a file.
 * 
 * @param n The first n coefficients of Fourier series.
 */
void Calculator::fourier(vector<double>::size_type n)
{
    // a is for Fourier's coefficient ak.
    // b is for Fourier's coefficient bk.
    // f is the function value at 2n+1 points.
    // t is the constant.
    // c and s are cos(2*pi/(2n+1)) and sin(2*pi/(2n+1)).
    // c1 and s1 are cos(kθ) and sin(kθ).
    // u1,u2,u0 are temporary variables.
    vector<double>::size_type i, j;
    vector<double> a, b, f;
    double t, c, s, c1, s1, u1, u2, u0;
    t = static_cast<double>(6.283185306 / (2.0 * n + 1.0));
    c = cos(t);
    s = sin(t);
    for (i = 0; i <= 2 * n; i++)
    {
        f.push_back(cal_num(orig, static_cast<double>((i + 0.5) * t)));
    }
    // Calculation process.
    t = static_cast<double>(2.0 / (2.0 * n + 1.0));
    c1 = 1.0;
    s1 = 0.0;
    for (i = 0; i <= n; i++)
    {
        u1 = 0.0;
        u2 = 0.0;
        for (j = 2 * n; j >= 1; j--)
        {
            u0 = f[j] + 2.0 * c1 * u1 - u2;
            u2 = u1;
            u1 = u0;
        }
        a.push_back(t * (f[0] + u1 * c1 - u2));
        b.push_back(t * u1 * s1);
        u0 = c * c1 - s * s1;
        s1 = c * s1 + s * c1;
        c1 = u0;
    }
    // Output operation.
    ofstream output;
    output.open("fourier_coefficients.txt");
    if (!output)
    {
        cout << "\nError: failed to open the input file." << endl;
        exit(0);
    }
    else
    {
        output << " i         "
               << "a           "
               << "   b" << endl;
        for (vector<double>::size_type i = 0; i < a.size(); i++)
        {
            output << setw(2) << i << "      " << setw(11) << a[i] << "       " << b[i] << endl;
        }
    }
    output.close();
    return;
}

/**
 * @brief Do Fast Fourier Transform.
 * 
 * @param k n=pow(2,k) is the number of sample points.
 * @param il il = 0 means not to calculate magnitude and angle; = 1 means to calculate them.
 * @param h Sample starting point.
 */
void Calculator::fft(vector<double>::size_type k, vector<double>::size_type il, double h)
{
    // pr is used to store the real part of input.
    // pi is used to store the imaginary part of input.
    // fr is used to return the real part of result.
    // fi is used to return the imaginary part of result.
    // it,m,is,i,j,nv,l0 are all for loop-iteration.
    // p,q,s,tr,ti,pdr,pdi are all for temporary purpose.
    vector<double>::size_type n = static_cast<vector<double>::size_type>(pow(2, k));
    vector<double> pr(n, 0), pi(n, 0), fr(n, 0), fi(n, 0), sq;
    vector<double>::size_type it, m, is, i, j, nv, l0;
    double p, q, s, tr, ti, pdr, pdi;

    // Get value of original points.
    for (i = 0; i < n; i++)
    {
        pr[i] = cal_num(orig, static_cast<double>(0.1 * i + 0.1 * h));
        pi[i] = 0.0;
    }
    sq = pr;

    // Get corresponding Even and odd parts.
    for (it = 0; it < n; it++)
    {
        m = it;
        is = 0;
        for (i = 0; i < k; i++)
        {
            j = m / 2;
            is = 2 * is + (m - 2 * j);
            m = j;
        }
        fr[it] = pr[is];
        fi[it] = pi[is];
    }

    // Do FFT to get the results.
    pr[0] = 1.0;
    pi[0] = 0.0;
    p = static_cast<double>(6.283185306 / (1.0 * n));
    pr[1] = cos(p);
    pi[1] = -sin(p);
    for (i = 2; i < n; i++)
    {
        p = pr[i - 1] * pr[1];
        q = pi[i - 1] * pi[1];
        s = (pr[i - 1] + pi[i - 1]) * (pr[1] + pi[1]);
        pr[i] = p - q;
        pi[i] = s - p - q;
    }
    for (it = 0; it < n - 1; it = it + 2)
    {
        tr = fr[it];
        ti = fi[it];
        fr[it] = tr + fr[it + 1];
        fi[it] = ti + fi[it + 1];
        fr[it + 1] = tr - fr[it + 1];
        fi[it + 1] = ti - fi[it + 1];
    }
    m = n / 2;
    nv = 2;
    for (l0 = k - 1; l0 >= 1; l0--)
    {
        m = m / 2;
        nv = 2 * nv;
        for (it = 0; it <= (m - 1) * nv; it = it + nv)
            for (j = 0; j <= (nv / 2) - 1; j++)
            {
                p = pr[m * j] * fr[it + j + nv / 2];
                q = pi[m * j] * fi[it + j + nv / 2];
                s = pr[m * j] + pi[m * j];
                s = s * (fr[it + j + nv / 2] + fi[it + j + nv / 2]);
                pdr = p - q;
                pdi = s - p - q;
                fr[it + j + nv / 2] = fr[it + j] - pdr;
                fi[it + j + nv / 2] = fi[it + j] - pdi;
                fr[it + j] = fr[it + j] + pdr;
                fi[it + j] = fi[it + j] + pdi;
            }
    }
    // Output operation.
    if (il == 1)
        for (i = 0; i < n; i++)
        {
            pr[i] = sqrt(fr[i] * fr[i] + fi[i] * fi[i]);
            if (fabs(fr[i]) < 0.000001 * fabs(fi[i]))
            {
                if ((fi[i] * fr[i]) > 0)
                    pi[i] = 90.0;
                else
                    pi[i] = -90.0;
            }
            else
                pi[i] = atan(fi[i] / fr[i]) * 360.0 / 6.283185306;
        }

    ofstream output;
    output.open("FFT_result.txt");
    if (!output)
    {
        cout << "\nError: failed to open the input file." << endl;
        exit(0);
    }
    else
    {
        output << " f(X)       " << setw(15) << "\"Real part after FFT\""
               << "   Imag part after FFT " << endl;
        for (vector<double>::size_type i = 0; i < n; i++)
        {
            output << setw(10) << sq[i] << setw(15) << fr[i] << "       " << setw(15) << fi[i] << endl;
        }
        output << " " << endl;

        if (il == 1)
        {
            output << " f(X)       " << setw(15) << "\"magnitude\""
                   << "               Angle" << endl;
            for (vector<double>::size_type i = 0; i < n; i++)
            {
                output << setw(10) << sq[i] << setw(15) << pr[i] << "       " << setw(15) << pi[i] << endl;
            }
        }
    }
    output.close();
    return;
}

/**
 * @brief Create random numbers.
 * 
 * @param u The average value.
 * @param g The variance of normal distribution. when g<=0, we ignore normal distribution restriction.
 * @param n The number of data.
 */
void Calculator::random(double u, double g, vector<double>::size_type n)
{
    // s, w, v are function coefficients listed in doc.
    // t and m for temporary usage.
    // r is the random seed.
    vector<double>::size_type i, k;
    vector<double> data(n, 0);
    int m;
    double s, w, v, t, r;
    time_t now = time(0);
    tm *ltm = localtime(&now);
    r = (ltm->tm_hour * ltm->tm_min) / ltm->tm_sec;
    s = 65536.0;
    w = 2053.0;
    v = 13849.0;
    if (g > 0)
    {
        for (k = 0; k < n; k++)
        {
            t = 0.0;
            for (i = 1; i <= 12; i++)
            {
                r = r * w + v;
                m = (int)(r / s);
                r = r - m * s;
                t = t + r / s;
            }
            data[k] = u + g * (t - 6.0);
        }
    }
    else
    {
        for (k = 0; k < n; k++)
        {
            m = (int)(r / s);
            r = r - m * s;
            r = w * r + v;
            m = (int)(r / s);
            r = r - m * s;
            data[k] = r / s;
        }
    }
    // Output operation
    ofstream output;
    output.open("random.txt");

    if (!output)
    {
        cout << "\nError: failed to open the input file." << endl;
        exit(0);
    }

    else
    {
        if (g > 0)
        {
            output << "Normal distribution random points with average= " << u << " variance= " << g << ":"
                   << "\n"
                   << endl;
            for (vector<double>::size_type i = 0; i < data.size(); i++)
            {
                output << setw(2) << data[i] << ",  ";
            }
            output << "End" << endl;
        }
        else
        {
            output << "random points in [0,1] are: "
                   << "\n"
                   << endl;
            for (vector<double>::size_type i = 0; i < data.size(); i++)
            {
                output << data[i] << ",  ";
            }
            output << "End" << endl;
        }
    }
    output.close();
    return;
}

/**
 * @brief Get the value of f(x).
 * 
 * @param a Lower boundary for solving value.
 * @param b Upper boundary for solving value.
 * @param h Step length.
 */
void Calculator::value(double a, double b, double h)
{
    double i = a;
    vector<double> t, ti;
    while (i < b)
    {
        t.push_back(cal_num(orig, i));
        ti.push_back(i);
        i = i + h;
    }
    for (vector<double>::size_type j = 0; j < t.size(); j++)
    {
        cout << "f(" << ti[j] << ")= " << t[j] << endl;
    }
}

/**
 * @brief Read input.txt and do corresponding calculation. Simplify the main file.
 * @param t Store temporary data.
 * @param str1 Store a whole row from input.txt.
 * @param Math_expr Store mathematical expression.
 * @param Deri_expr Store derivative expression.
 * @param res Arithmetic calculation result.
 * @param Opr_type Supported calculations.
 */
void cal_expression()
{
    string t;
    string str1;
    string Math_expr;
    string Deri_expr = "0";
    basic_string<char>::size_type l, a;
    double res;
    vector<double> parameter, zeros;
    enum Opr_type
    {
        Limit = 1,
        Integration,
        Fourier_coefficients,
        Get_zero,
        Integration_0_inf,
        FFT,
        Random,
        Value
    };

    map<string, Opr_type> Opr_str{{"Limit", Limit}, {"Integration", Integration}, {"Fourier_coefficients", Fourier_coefficients}, {"Get_zero", Get_zero}, {"Integration_0_inf", Integration_0_inf}, {"FFT", FFT}, {"Random", Random}, {"Value", Value}};

    // Read mathematical expression from input.txt and detect whether the file is open or not.
    ifstream input("input.txt");
    if (!input.is_open())
    {
        cout << "\nError: failed to open the input file." << endl;
        exit(1);
    }

    // First loop to read the expression.
    while (getline(input, str1))
    {

        // Remove spaces in str
        str1.erase(remove(str1.begin(), str1.end(), ' '), str1.end());

        // Extract the original expression from file into Math_expre.
        if (str1.find("f(X)") < str1.length() && str1.find('#') == string::npos)
        {
            Math_expr = str1.substr(str1.find("f(X)") + 5);
            cout << "f(X)=" << Math_expr << endl;
            l = Math_expr.length();
            a = 0;
            while (a < l)
            {
                if (Math_expr[a] == '-' && (a == 0 || Math_expr[a - 1] == '('))
                {
                    Math_expr.replace(a, 1, "0-");
                    a = a + 2;
                    l = l + 1;
                }
                else
                {
                    a = a + 1;
                    continue;
                }
            }
            // If there is no variable in expression, then just do arithmetic calculation.
            if (Math_expr.find('X') == string::npos)
            {
                Calculator cal(Math_expr, Deri_expr);
                cal.getVal(res);
                cout << "f(X)" << "=" << res << endl;
                return;
            }
        }
        // Extract the derivative expression of expr.
        if (str1.find("f'(X)") < str1.length() && str1.find('#') == string::npos)
        {
            Deri_expr = str1.substr(str1.find("f'(X)") + 6);
            cout << "f'(X)=" << Deri_expr << "\n"
                 << endl;
            l = Deri_expr.length();
            a = 0;
            while (a < l)
            {
                if (Deri_expr[a] == '-' && (a == 0 || Deri_expr[a - 1] == '('))
                {
                    Deri_expr.replace(a, 1, "0-");
                    a = a + 2;
                    l = l + 1;
                }
                else
                {
                    a = a + 1;
                    continue;
                }
            }
        }
        if (str1.find("Operations") < str1.length())
            break;
    }

    Calculator cal(Math_expr, Deri_expr);

    // Second loop to do specific mathematical operation.
    while (getline(input, str1))
    {
        if (str1.find("End") < str1.length())
            break;

        if (!(str1.empty()) && str1.find('#') == string::npos)
        {
            stringstream tem(str1);
            tem >> t;
            Opr_type type = Opr_str[t];
            if (type == 0)
            {
                cout << "\nError: invalid operation, please check your input.\n"
                     << endl;
                exit(2);
            }
            switch (type)
            {
            case Limit:
                cout << "***** Limit Operation *****" << endl;
                while (tem >> t)
                {
                    if (!cal.num_detect(t))
                    {
                        cout << "\nError: invalid input format for limit. It should be a number." << endl;
                        exit(2);
                    }
                    parameter.push_back(stod(t));
                }
                if (parameter.size() != 4)
                {
                    cout << "*****Error: too many or not enough data for limit.*****" << endl;
                    exit(2);
                }
                if (parameter[0]>=parameter[1])
                {
                    cout << "Error: upper boundary is smaller than lower boundary." << endl;
                    exit(2);
                }
                if (parameter[2]<0.001 || parameter[2]>abs(parameter[0]-parameter[1]))
                {
                    cout << "Error: Step length for limit is unsuitable." << endl;
                    exit(2);
                }
                cal.limit(parameter[0], parameter[1], parameter[2], parameter[3]);
                parameter.clear();
                break;

            case Integration:
                cout << "***** Integration Operation *****" << endl;
                while (tem >> t)
                {
                    if (!cal.num_detect(t))
                    {
                        cout << "\nError: invalid input format for integration. It should be a number." << endl;
                        exit(2);
                    }
                    parameter.push_back(stod(t));
                }
                if (parameter.size() != 3)
                {
                    cout << "Error: too many or not enough data for integration." << endl;
                    exit(2);
                }
                if (parameter[0]>=parameter[1])
                {
                    cout << "Error: upper boundary is smaller than lower boundary." << endl;
                    exit(2);
                }
                parameter.push_back(cal.integra(Math_expr, parameter[0], parameter[1], parameter[2]));
                cout << "integration from " << parameter[0] << " to " << parameter[1] << " of f(X)"
                     << " is: " << parameter[3] << "\n"
                     << endl;
                parameter.clear();
                break;

            case Fourier_coefficients:
                cout << "***** Fourier_coefficients Operation *****" << endl;
                while (tem >> t)
                {
                    if (!cal.num_detect(t))
                    {
                        cout << "\nError: invalid input format for Fourier_transform. It should be a number." << endl;
                        exit(2);
                    }
                    parameter.push_back(stod(t));
                }
                if (parameter.size() != 1)
                {
                    cout << "Error: too many or not enough data for Fourier_transform." << endl;
                    exit(2);
                }
                if (parameter[0] < 0)
                {
                    cout << "Error: input number for Fourier_transform should be a positive number." << endl;
                    exit(2);
                }
                cal.fourier(static_cast<vector<double>::size_type>(parameter[0]));
                cout << "Successfully calculate Fourier_coefficients and store results in fourier_coefficients.txt.\n";
                parameter.clear();
                break;

            case Get_zero:
                cout << "***** Get zero Operation *****" << endl;
                while (tem >> t)
                {
                    if (!cal.num_detect(t))
                    {
                        cout << "\nError: invalid input format for Get zero. It should be a number." << endl;
                        exit(2);
                    }
                    parameter.push_back(stod(t));
                }
                if (parameter.size() != 4)
                {
                    cout << "Error: too many or not enough data for Get zero." << endl;
                    exit(2);
                }
                if (parameter[0]>=parameter[1])
                {
                    cout << "Error: upper boundary is smaller than lower boundary." << endl;
                    exit(2);
                }
                if (parameter[2]<0.001 || parameter[2]>abs(parameter[0]-parameter[1]))
                {
                    cout << "Error: Step length for zero_operation is unsuitable." << endl;
                    exit(2);
                }
                zeros = cal.zero(Math_expr, parameter[0], parameter[1], parameter[2], parameter[3]);
                parameter.clear();
                break;

            case Integration_0_inf:
                if (Math_expr.find("exp(-X)") == string::npos && Math_expr.find("exp(0-X)") == string::npos)
                {
                    cout << "Error: the form of Integration_0_inf should be exp(-X)*f(X)!" <<endl;
                    exit(3);
                }
                cout << "***** Integrate from 0 to infinity *****" << endl;
                while (tem >> t)
                {
                    cout << "Error: too many arguments for integration from 0 to infinity." << endl;
                    exit(2);
                }
                cout << "integration from 0"
                     << " to "
                     << " infinity of f(X)"
                     << " is: " << cal.integra_0i(Math_expr) << "\n"
                     << endl;
                break;

            case FFT:
                cout << "\n***** FFT Operation *****" << endl;
                while (tem >> t)
                {
                    if (!cal.num_detect(t))
                    {
                        cout << "\nError: invalid input format for Get zero. It should be a number." << endl;
                        exit(2);
                    }
                    parameter.push_back(stod(t));
                }
                if (parameter.size() != 3)
                {
                    cout << "Error: too many or not enough data for Get zero." << endl;
                    exit(2);
                }
                if (parameter[0] <= 1 || !(parameter[1] == 0 || parameter[1] == 1))
                {
                    cout << "error: input for FFT have unsuitable value." << endl;
                    exit(3);
                }
                cal.fft(static_cast<vector<double>::size_type>(parameter[0]), static_cast<vector<double>::size_type>(parameter[1]), parameter[2]);
                cout << "Successfully calculate FFT and store results in FFT_result.txt.\n";
                parameter.clear();
                break;

            case Random:
                cout << "\n***** random Operation *****" << endl;
                while (tem >> t)
                {
                    if (!cal.num_detect(t))
                    {
                        cout << "\nError: invalid input format for random. It should be a number." << endl;
                        exit(2);
                    }
                    parameter.push_back(stod(t));
                }

                if (parameter.size() != 3)
                {
                    cout << "Error: too many or not enough data for Get zero." << endl;
                    exit(2);
                }

                if (parameter[1]<0 || parameter[2] < 10)
                {
                    cout << "Error: invalid input parameter for random operation." << endl;
                    exit(3);
                }
                cal.random(parameter[0], parameter[1], static_cast<vector<double>::size_type>(parameter[2]));
                cout << "Successfully calculate random points and store results in random.txt.\n";
                parameter.clear();
                break;

            case Value:
                cout << "\n***** value Operation *****" << endl;
                while (tem >> t)
                {
                    if (!cal.num_detect(t))
                    {
                        cout << "\nError: invalid input format for calculating value. It should be a number." << endl;
                        exit(2);
                    }
                    parameter.push_back(stod(t));
                }

                if (parameter.size() != 3)
                {
                    cout << "Error: too many or not enough data for value." << endl;
                    exit(2);
                }

                if (parameter[0]>=parameter[1])
                {
                    cout << "Error: upper boundary is smaller than lower boundary." << endl;
                    exit(2);
                }
                if (parameter[2]<0.001 || parameter[2]>abs(parameter[0]-parameter[1]))
                {
                    cout << "Error: Step length for getting value is unsuitable." << endl;
                    exit(2);
                }
                cal.value(parameter[0], parameter[1], parameter[2]);
                parameter.clear();
                break;
            }
        }
    }
    input.close();
}
