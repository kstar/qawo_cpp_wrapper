/***************************************************************************
                         qawo_cpp_wrapper.cpp  -
                             -------------------
    begin                : Wed 19 Oct 2016 17:56:04 CDT
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


/* Project Includes */
#include "qawo_cpp_wrapper.h"
#include <assert.h>

double callIntegrand( double x, void *params ) {
    return reinterpret_cast<Quadrature *>( params )->m_integrand( x );
}

Quadrature::Quadrature() {
    F.function = &callIntegrand;
    F.params = this;
    m_limit = 1000;
    m_bisections = 24;
    m_epsabs = 1e-7;
    m_epsrel = 1e-4;
}

std::vector<double> Quadrature::transform( std::function<double( double )> integrand , const std::vector<double>& omegaVector, double a, double b, const gsl_integration_qawo_enum trigType ) {
    std::vector<double> result;
    bool infiniteInterval;
    double L;

    m_errors.clear();

    if ( omegaVector.empty() )
        return result; // nothing to do

    if ( a == b )
        return result; // nothing to do. FIXME: Technically we should return a vector of zeros

    /* Figure out whether we have a finite or semi-infinite interval */
    if ( !std::isfinite( a ) || !std::isfinite( b ) )
        infiniteInterval = true;

    if ( !std::isfinite( a ) ) {
        // FIXME: Handle this case
        std::cerr << "Error: This method can only compute semi-infinite intervals of the form [a, +inf) as of this writing. Manually transform your integral to the required form." << std::endl;
        return result;
    }

    if ( b <= a ) {
        // FIXME: Handle this case, maybe?
        std::cerr << "Error: upper limit of integral must be larger than lower limit as of this writing. Manually flip your integral to the required form." << std::endl;
        return result;
    }

    /* Set the integrand */
    m_integrand = std::move( integrand );

    /* Allocate memory for the results */
    m_errors.reserve( omegaVector.size() );
    result.reserve( omegaVector.size() );

    /* Allocate workspace */
   gsl_integration_workspace *workspace = gsl_integration_workspace_alloc( m_limit );
   gsl_integration_workspace *cycle_workspace = NULL;

    /* Allocate the QAWO table */
    L = ( infiniteInterval ? 1. : b - a );
    gsl_integration_qawo_table *table = gsl_integration_qawo_table_alloc( omegaVector[0], L, trigType, m_bisections );

    if ( infiniteInterval )
        cycle_workspace = gsl_integration_workspace_alloc( m_limit );

    for ( double omega : omegaVector ) {
        double singleResult, abserr;
        gsl_integration_qawo_table_set( table, omega, L, trigType );
        if ( infiniteInterval )
            gsl_integration_qawf( &F, a, m_epsabs, m_limit, workspace, cycle_workspace, table, &singleResult, &abserr );
        else
            gsl_integration_qawo( &F, a, m_epsabs, m_epsrel, m_limit, workspace, table, &singleResult, &abserr );
        result.push_back( singleResult );
        m_errors.push_back( abserr );
    }

    /* Free memory */
    gsl_integration_qawo_table_free( table );
    gsl_integration_workspace_free( workspace );
    if ( infiniteInterval )
        gsl_integration_workspace_free( cycle_workspace );

    return result;
}
