#pragma once

inline PS::F64 mutualRadius(PS::F64 r1, PS::F64 r2){
    return std::max(r1, r2);
}

inline PS::F64 cutoff_f(PS::F64 y, PS::F64 g, PS::F64 ginv){
    //const PS::F64 g = FPGrav::gamma;
    //const PS::F64 ginv = pow(g-1., 7.);
    
    PS::F64 y2 = y*y;
    PS::F64 y3 = y2*y;
    PS::F64 y4 = y2*y2;
    PS::F64 y5 = y3*y2;
    PS::F64 y6 = y3*y3;
    PS::F64 y7 = y4*y3;
    
    PS::F64 f = 0.0;
    f =( -10.0*y7/3.0 +14.0*(g+1.0)*y6 -21.0*(g*g +3.0*g +1.0)*y5
         +y4*35.0*(g*g*g +9.0*g*g +9.0*g +1.0)/3.0 -70.0*y3*(g*g*g +3.0*g*g + g)
         +210.0*(g +1.0)*g*g*y2 -140.0*g*g*g*y*log(y) +g*g*g*g*(g*g*g -7.0*g*g +21.0*g -35.0) )
        / ginv;

    return f;   
}

//#ifndef FORDEBUG

inline PS::F64 cutoff_W(PS::F64 rij, PS::F64 r_out){
    const PS::F64 g = FPGrav::gamma;
    static const PS::F64 ginv = pow(FPGrav::gamma-1., 7.);
    PS::F64 y = rij/r_out;

    PS::F64 W = 0.0;
    if ( 1.0 <= y ) {
        W = 1.0;
    } else if ( y <= g ){
        PS::F64 g2 = g*g;
        PS::F64 g4 = g2*g2;
        PS::F64 g5 = g4*g;
        PS::F64 g6 = g4*g2;
        W = 7.0*y*(g6-9.*g5+45.*g4-60.*g2*g*log(g)-45.*g2+9.*g-1.)/(3.*ginv);
    } else {
        PS::F64 f1 = cutoff_f(1., g, ginv);
        PS::F64 f = cutoff_f(y, g, ginv);
        W = f +y*(1.-f1);
    }

    return W;
}

inline PS::F64 cutoff_K(PS::F64 rij, PS::F64 r_out){
    const PS::F64 g = FPGrav::gamma;
    PS::F64 y = rij/r_out;
    PS::F64 x = (y - g)/(1. - g);
    
    PS::F64 K = 0.;
    if( x < 0. ) {
        K = 0.;
    } else if ( x >= 1.0 ) {
        K = 1.;
    } else {
        PS::F64 x2 = x*x;
        PS::F64 x4 = x2*x2;
        K = (((-20.0*x+70.0)*x-84.0)*x+35.0)*x4;
    }
    return K;
}

inline PS::F64 cutoff_dK(PS::F64 rij, PS::F64 rijvij, PS::F64 r_out){
    const PS::F64 g = FPGrav::gamma;
    PS::F64 y = rij/r_out;
    PS::F64 x = (y - g)/(1. - g);
    PS::F64 dx = rijvij/(r_out*rij*(1.-g));

    PS::F64 dK = 0.;
    if( x < 0. || x >= 1.0 ) {
        dK = 0.;
    } else {
        PS::F64 x2 = x*x;
        PS::F64 x3 = x2*x;
        PS::F64 x4 = x2*x2;
        PS::F64 x5 = x4*x;
        PS::F64 x6 = x4*x2;
        dK = (-140.0*x6 + 420.0*x5 - 420.0*x4 + 140.0*x3) * dx;
    }
    return dK;
}

//#else
//inline PS::F64 cutoff_W(PS::F64 rij, PS::F64 r_out){ return 0.; }
//inline PS::F64 cutoff_K(PS::F64 rij, PS::F64 r_out){ return 0.; }
//inline PS::F64 cutoff_dK(PS::F64 rij, PS::F64 rijvij, PS::F64 r_out){ return 0.;}
//#endif

inline PS::F64 cutoff_W2(PS::F64 dr2, PS::F64 r_out){
    return cutoff_W(sqrt(dr2), r_out);
}

inline PS::F64 cutoff_W2(PS::F64 dr2, PS::F64 r_out1, PS::F64 r_out2){
    return cutoff_W(sqrt(dr2), mutualRadius(r_out1, r_out2));
}
