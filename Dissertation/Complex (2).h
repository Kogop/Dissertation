// Complex number library

#ifndef __cplusplus
#error Must use C++ for the type complex.
#endif

#ifndef _COMPLEX_H
#define _COMPLEX_H
//#include "data.h"
#if !defined( __MATH_H )
#include <math.h>
/* Defined for compatibility with foolish Microsoft math.h
  Definition:
  #define complex _complex
  presents in math.h. It's really fool, isn't it?
*/

#undef complex
#define complex complex
#endif

class complex {
public:
	// constructors
	complex(double __re_val, double __im_val = 0);
	complex();

	// complex manipulations
	friend double real(const complex&);  // the real part
	friend double imag(const complex&);  // the imaginary part
	friend complex conj(const complex&); // the complex conjugate
	friend double norm(const complex&);  // the square of the magnitude
	friend double arg(const complex&);  // the angle in the plane
	friend complex Root(const complex&); // square root

	// Overloaded ANSI C math functions
	friend double  abs(const complex&);
	friend complex cos(const complex&);
	friend complex exp(const complex&);
	friend complex sin(const complex&);
	

	// Binary Operator Functions
	friend complex operator+(const complex&, const complex&);
	friend complex operator+(double, const complex&);
	friend complex operator+(const complex&, double);
	friend complex operator-(const complex&, const complex&);
	friend complex operator-(double, const complex&);
	friend complex operator-(const complex&, double);
	friend complex operator*(const complex&, const complex&);
	friend complex operator*(const complex&, double);
	friend complex operator*(double, const complex&);
	friend complex operator/(const complex&, const complex&);
	friend complex operator/(const complex&, double);
	friend complex operator/(double, const complex&);
	friend int operator==(const complex&, const complex&);
	friend int operator!=(const complex&, const complex&);
	const complex& operator+=(const complex&);
	const complex& operator+=(double);
	const complex& operator-=(const complex&);
	const complex& operator-=(double);
	const complex& operator*=(const complex&);
	const complex& operator*=(double);
	const complex& operator/=(const complex&);
	const complex& operator/=(double);
	complex operator+();
	complex operator-();

	// Implementation
private:
	double re, im;
};

// i constant
const complex _i = complex(0, 1);

// Inline complex functions

inline complex::complex(double __re_val, double __im_val)
{
	re = __re_val;
	im = __im_val;
}

inline complex::complex()
{
	/* if you want your complex numbers initialized ...
	   Yes, I want it! */
	re = im = 0;
}

inline complex complex::operator+()
{
	return *this;
}

inline complex complex::operator-()
{
	return complex(-re, -im);
}

// Definitions of compound-assignment operator member functions

inline const complex& complex::operator+=(const complex& __z2)
{
	re += __z2.re;
	im += __z2.im;
	return *this;
}

inline const complex& complex::operator+=(double __re_val2)
{
	re += __re_val2;
	return *this;
}

inline const complex& complex::operator-=(const complex& __z2)
{
	re -= __z2.re;
	im -= __z2.im;
	return *this;
}

inline const complex& complex::operator-=(double __re_val2)
{
	re -= __re_val2;
	return *this;
}

inline const complex& complex::operator*=(const complex& __z2)
{
	double __re_val1 = re;
	re = __re_val1 * __z2.re - im * __z2.im;
	im = __re_val1 * __z2.im + im * __z2.re;
	return *this;
}

inline const complex& complex::operator*=(double __re_val2)
{
	re *= __re_val2;
	im *= __re_val2;
	return *this;
}

inline const complex& complex::operator/=(double __re_val2)
{
	re /= __re_val2;
	im /= __re_val2;
	return *this;
}

inline const complex& complex::operator/=(const complex& __z2)
{
	double __norm = __z2.re * __z2.re + __z2.im * __z2.im;
	double __re_val1 = re;
	re = (__re_val1 * __z2.re + im * __z2.im) / __norm;
	im = (im * __z2.re - __re_val1 * __z2.im) / __norm;
	return *this;
}

// Definitions of non-member complex functions

inline double real(const complex& __z)
{
	return __z.re;
}

inline double imag(const complex& __z)
{
	return __z.im;
}

inline complex conj(const complex& __z)
{
	return complex(__z.re, -__z.im);
}

// Definitions of non-member binary operator functions

inline complex operator+(const complex& __z1, const complex& __z2)
{
	return complex(__z1.re + __z2.re, __z1.im + __z2.im);
}

inline complex operator+(double __re_val1, const complex& __z2)
{
	return complex(__re_val1 + __z2.re, __z2.im);
}

inline complex operator+(const complex& __z1, double __re_val2)
{
	return complex(__z1.re + __re_val2, __z1.im);
}

inline complex operator-(const complex& __z1, const complex& __z2)
{
	return complex(__z1.re - __z2.re, __z1.im - __z2.im);
}

inline complex operator-(double __re_val1, const complex& __z2)
{
	return complex(__re_val1 - __z2.re, -__z2.im);
}

inline complex operator-(const complex& __z1, double __re_val2)
{
	return complex(__z1.re - __re_val2, __z1.im);
}

inline complex operator*(const complex& __z1, double __re_val2)
{
	return complex(__z1.re * __re_val2, __z1.im * __re_val2);
}

inline complex operator*(double __re_val1, const complex& __z2)
{
	return complex(__z2.re * __re_val1, __z2.im * __re_val1);
}

inline complex operator*(const complex& __z1, const complex& __z2)
{
	return complex(__z1.re * __z2.re - __z1.im * __z2.im,
		__z1.re * __z2.im + __z1.im * __z2.re);
}

inline complex operator/(const complex& __z1, double __re_val2)
{
	return complex(__z1.re / __re_val2, __z1.im / __re_val2);
}

inline complex operator/(const complex& __z1, const complex& __z2)
{
	double __norm = __z2.re * __z2.re + __z2.im * __z2.im;
	return complex((__z1.re * __z2.re + __z1.im * __z2.im) / __norm,
		(__z1.im * __z2.re - __z1.re * __z2.im) / __norm);
}

inline complex operator/(double __re_val_1, const complex& __z2)
{
	double __norm = __z2.re * __z2.re + __z2.im * __z2.im;
	return complex(__re_val_1 * __z2.re / __norm,
		-__re_val_1 * __z2.im / __norm);
}

inline int operator==(const complex& __z1, const complex& __z2)
{
	return __z1.re == __z2.re && __z1.im == __z2.im;
}

inline int operator!=(const complex& __z1, const complex& __z2)
{
	return __z1.re != __z2.re || __z1.im != __z2.im;
}

// Overloaded ANSI C math functions
inline double abs(const complex& __z)
{
	return sqrt(__z.re * __z.re + __z.im * __z.im);
}

inline complex cos(const complex& __z)
{
	double __e_im1 = exp(__z.im) / 2,
		__e_im2 = exp(-__z.im) / 2;
	return complex(cos(__z.re) * (__e_im1 + __e_im2),
		sin(__z.re) * (__e_im2 - __e_im1));
}

inline complex sin(const complex& __z)
{
	double __e_im1 = exp(__z.im) / 2,
		__e_im2 = exp(-__z.im) / 2;
	return complex(sin(__z.re) * (-__e_im1 - __e_im2),
		cos(__z.re) * (__e_im1 - __e_im2));
}

inline complex exp(const complex& __z)
{
	double __exp = exp(__z.re);
	return complex(__exp * cos(__z.im), __exp * sin(__z.im));
}

//complex complex::Root(complex z) {
//	double fi;
//	complex a1, a2;
//	fi = 2 * atan((z.imag) / (z.real + abs(z)));
//	a1 = (sqrt(abs(z)) * (cos((fi) / 2) + _i * sin((fi) / 2)));
//	return a1;
//}

// Complex stream I/O

//#include <fstream>

//ostream& operator<<(ostream&, const complex&);
//istream& operator>>(istream&, complex&);

#endif // __COMPLEX_H
