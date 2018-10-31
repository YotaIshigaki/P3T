#pragma once


PS::F64 pow2(PS::S32 x)
{
    if ( x == 0 ){
        return 1.;
    } else if ( x > 0 ) {
        return 2.*pow2(x-1);
    } else {
        return 0.5*pow2(x+1);
    }
}

PS::F64 getvalue(std::string value,
                 PS::F64 MKS_UNIT,
                 PS::F64 CGS_UNIT)
{
    PS::F64 result = 1.;
    std::vector<std::string> numer, denom;
    numer.clear();
    denom.clear();

    if ( value.size() > 4 ){
        std::string last3 = value.substr(value.size()-3, 3);
        if ( last3 == "MKS" || last3 == "mks" ) {
            result /= MKS_UNIT;
            value.erase(value.end()-3, value.end());
        }
        if ( last3 == "CGS" || last3 == "cgs" ) {
            result /= CGS_UNIT;
            value.erase(value.end()-3, value.end());
        }
    }

    PS::S32 begin = 0;
    PS::F64 x = 0;
    bool isNumer = true;
    for ( PS::S32 i=0; i<value.size()+1; i++ ){
        if ( i == value.size() || value[i] == '*' || value[i] == '/' ){
            std::string x_str = value.substr(begin, i-begin);
            
            if ( x_str.size() > 2 && x_str.substr(0, 2) == "2^" ) {
                x = pow2( std::atoi( x_str.substr(2, x_str.size()-2).c_str() ) );
            } else {
                x = std::atof( x_str.c_str() );
            }
            begin = i+1;

            if ( isNumer ) {
                result *= x;
            } else {
                result /= x;
            }
            
            if ( value[i] == '*' ) isNumer = true;
            if ( value[i] == '/' ) isNumer = false;
        }
    }
    
    return result;
}

PS::F64 KeplerEq(PS::F64 u,
                 PS::F64 ecc){ return u - ecc * sin(u); }
PS::F64 solveKeplerEq(PS::F64 l,
                      PS::F64 ecc)
{
    PS::F64 u;
    
    PS::F64 ecc2 = ecc*ecc;
    PS::F64 ecc3 = ecc2*ecc;
    PS::F64 ecc4 = ecc2*ecc2;
    PS::F64 ecc5 = ecc3*ecc2;
    PS::F64 ecc6 = ecc3*ecc3;
    u = l
        + (ecc     - ecc3/8. + ecc5/192.)*sin(l)
        + (ecc2/2. - ecc4/6. + ecc6/48. )*sin(2.*l)
        + (3.*ecc3/8. - 27.*ecc5/128.)*sin(3.*l)
        + (   ecc4/3. -  4.*ecc6/15. )*sin(4.*l)
        + 125.*ecc5/384.*sin(5.*l)
        +  27.*ecc6/ 80.*sin(6.*l);
#if 1
    if ( abs(KeplerEq(u,ecc)-l) > 1.e-15  ){
        PS::F64 u0;
        PS::S32 loop = 0;
        //u = l;
        do {
            u0 = u;
            PS::F64 sinu0 = sin(u0);
            PS::F64 cosu0 = cos(u0);
            u = u0 - ((u0 - ecc*sinu0 - l)/(1. - ecc*cosu0));
            loop++;
        } while( fabs(u - u0) > 1.e-15 && loop < 10 );
    }
#endif

    return u;
}
void posVel2OrbitalElement(PS::F64vec pos,
                           PS::F64vec vel,
                           PS::F64 mu,
                           PS::F64 & ax,
                           PS::F64 & ecc,
                           PS::F64 & n,
                           PS::F64 & u,
                           PS::F64vec & P,
                           PS::F64vec & Q)
{
    PS::F64 r2 = pos*pos;
    PS::F64 r = sqrt(r2);
    PS::F64 rinv = 1./r;
    PS::F64 v2 = vel*vel;
    PS::F64 rv = pos*vel;
    
    ax = 1.0 / (2.0*rinv - v2 / mu);
    PS::F64 ecccosu = 1. - r/ax;
    PS::F64 eccsinu = rv/sqrt(mu*ax);
    ecc = sqrt(ecccosu*ecccosu + eccsinu*eccsinu);
    n = sqrt(mu / (ax*ax*ax));
    
    u = atan2(eccsinu, ecccosu);
    
    PS::F64 cosu = ecccosu/ecc;
    PS::F64 sinu = eccsinu/ecc;
    PS::F64 aninv = sqrt(ax/mu);
    PS::F64 ecc_sq = sqrt(1.-ecc*ecc);
    P = rinv*cosu * pos - aninv*sinu * vel;
    Q = (rinv*sinu * pos + aninv*(cosu-ecc) * vel )/ecc_sq;
}
void orbitalElement2PosVel(PS::F64vec & pos,
                           PS::F64vec & vel,
                           PS::F64 mu,
                           PS::F64 ax,
                           PS::F64 ecc,
                           PS::F64 n,
                           PS::F64 u,
                           PS::F64vec P,
                           PS::F64vec Q)
{
    PS::F64 cosu = cos(u);
    PS::F64 sinu = sin(u);
    PS::F64 ecc_sq = sqrt(1.-ecc*ecc);
    pos = ax * ((cosu-ecc) * P + ecc_sq*sinu * Q);
    PS::F64 rinv = sqrt(1./(pos*pos));
    vel = ax*ax*n*rinv* (-sinu * P + ecc_sq*cosu * Q);
}
