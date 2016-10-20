/***************************************************************************
                         qawo_wrapper_test.cpp  -
                             -------------------
    begin                : Wed 19 Oct 2016 19:39:28 CDT
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

/* STL Includes */
#include <cmath>
#include <assert.h>

int main() {
    std::vector<double> t;
    for ( double p = -3.; p <= +3.; p += 0.1 )
        t.push_back( std::pow( 10, p ) );
    std::cout << "# t, f(t), (pi/2) * exp(-t)" << std::endl;
    Quadrature q;
    /* Transform a Lorentzian and check it against a decaying exponential */
    std::vector<double> f = q.cosineTransform( []( double w ) -> double { return 1./( std::pow( w, 2 ) + 1 ); }, t );
    assert( f.size() == t.size() );
    for ( int i = 0; i < t.size(); ++i )
        std::cout << t[i] << ", " << f[i] << ", " << ( M_PI/2. ) * std::exp( -t[i] ) << std::endl;
}
