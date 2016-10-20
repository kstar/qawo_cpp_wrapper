/***************************************************************************
                          qawo_cpp_wrapper.h  -  
                             -------------------
    begin                : Wed 19 Oct 2016 17:39:21 CDT
    copyright            : (c) 2016 by Akarsh Simha
    email                : akarsh.simha@kdemail.net
***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/



#ifndef QAWO_CPP_WRAPPER_H
#define QAWO_CPP_WRAPPER_H

#include <gsl/gsl_integration.h>

#include <functional>
#include <vector>
#include <limits>
#include <iostream>


/**
 * @class Quadrature
 * @short A C++ wrapper for GSL's routines for quadrature against cos(kx) or sin(kx). (i.e. oscillatory integrals)
 * @author Akarsh Simha <akarsh.simha@kdemail.net>
 */

class Quadrature {

public:

    /**
     * @short Constructor
     */
    Quadrature();

    /**
     * @short Sets the relative and absolute error bounds desired
     */
    inline void setErrorBounds( double epsabs, double epsrel ) { m_epsabs = epsabs; m_epsrel = epsrel; }

    /**
     * @short Sets the limiting number of bisections for the Tchebyshev moment table
     */
    inline void setQawoTableBisections( std::size_t N ) { m_bisections = N; }

    /**
     * @short Set subinterval limit
     * Sets the limit parameter of the QAWO / QAWF routines; and also uses the same limit for table allocation
     */
    inline void setSubintervalLimit( const std::size_t limit ) { m_limit = limit; }

    /**
     * @short Sine integral
     * Computes \int_{a}^{b} integrand(x) sin( omega x ) dx by calling GSL routines QAWO and QAWF
     * @note b can be inf
     */
    inline std::vector<double> sineTransform( std::function<double(double)> integrand, const std::vector<double> &omegaVector, const double a = 0, const double b = std::numeric_limits<double>::infinity() ) { return transform( std::move(integrand), omegaVector, a, b, GSL_INTEG_SINE ); }

    /**
     * @short Cosine integral
     * @sa sineTransform for further details
     */
    inline std::vector<double> cosineTransform( std::function<double(double)> integrand, const std::vector<double> &omegaVector, const double a = 0, const double b = std::numeric_limits<double>::infinity() ) { return transform( std::move(integrand), omegaVector, a, b, GSL_INTEG_COSINE ); }

    /**
     * @short Returns the vector of absolute errors for each value of omega, as given returned by QAWO / QAWF
     */
    inline const std::vector<double> &lastIntegralErrorVector() const { return m_errors; }

private:

    friend double callIntegrand( double x, void *params );

    /**
     * @short The function that does all the actual work
     */
    std::vector<double> transform( std::function<double(double)> integrand, const std::vector<double> &omegaVector, const double a, const double b, const enum gsl_integration_qawo_enum trigType );

    /**
     * @short Set-up for integration
     */
    std::function<double(double)> m_integrand;
    std::vector<double> m_errors;
    double m_epsabs, m_epsrel;
    gsl_function F;
    std::size_t m_limit, m_bisections;
};

/**
 * @short A wrapper to call a std::function object, so that we can
 * get a function pointer to it, as GSL requires that.
 */
#endif
