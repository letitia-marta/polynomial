# Polynomial Calculator with Multiple Coefficient Types ✨

This project is a **C++ console application** for creating and manipulating **polynomials** with different types of coefficients:  
- **Integers**  
- **Real numbers (float)**  
- **Rational numbers**  
- **Complex numbers**

It supports polynomial operations (addition, multiplication, evaluation, derivative) and implements **Newton's Method** for root finding in the real case.

## 🚀 Features

### 🔹 Supported Coefficient Types
- `int` → Integer coefficients  
- `float` → Real number coefficients (includes Newton’s root finding)  
- `NrRational` → Custom rational numbers class  
- `NrComplex` → Custom complex numbers class  

### 🔹 Polynomial Operations
- **Addition (`+`)**  
- **Multiplication (`*`)**  
- **Evaluation (`[]`)** – compute polynomial value at a given point  
- **Derivative** computation  

### 🔹 Extra Features
- **Newton’s Method** for root approximation on real polynomials:  xn+1 = xn - f(xn)/f'(xn)
- Overloaded operators for rational and complex arithmetic.  
- Clear and structured output formatting.  

## 🖥️ Usage

When running the program, a menu will appear:
- Choose a type of coefficient (1–4).  
- Enter the degree and coefficients for two polynomials.  
- The program displays their **sum**, **product**, and evaluates them for a given input.  
- For real polynomials, you can also provide an interval `[a, b]` to search for roots using **Newton’s method**.  

## 📊 Example Inputs

### Integers
3

2 3 0 -1

2

3 2 1

Polynomials:
- \( f(x) = 2x^{3} + 3x^{2} - 1 \)  
- \( g(x) = 3x^{2} + 2x + 1 \)  

### Reals
2

1 -2 1

3

2.1 3 0 -1.4

Polynomials:
- \( f(x) = x^{2} - 2x + 1 \)  
- \( g(x) = 2.1x^{3} + 3x^{2} - 1.4 \)  

### Rationals
2

2 1 3 2 1 1

1

3 1 2 4

Polynomials:
- \( f(x) = 2x^{2} + \tfrac{3}{2}x + 1 \)  
- \( g(x) = 3x + \tfrac{2}{4} \)  

### Complex
2

2 1 3 2 1 1

1

3 1 2 4

Polynomials:
- \( f(x) = (2+i)x^{2} + (3+2i)x + (1+i) \)  
- \( g(x) = (3+i)x + (2+4i) \)  

## ⚙️ Implementation Details
- **Classes**:
  - `NrRational` → Rational number arithmetic with reduction using `cmmdc` (greatest common divisor).  
  - `NrComplex` → Complex number arithmetic with custom operators.  
  - `Polinom<T>` → Templated polynomial class supporting multiple coefficient types.  

- **Key Functions**:
  - `alg<T>(f, g)` → Reads two polynomials, outputs sum, product, evaluations.  
  - `Newton(f, a, b)` → Applies Newton’s method for root approximation.  
  - `clearscreen()` → Clears console between runs.
