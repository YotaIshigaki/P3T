#pragma once

///////////////////////////////////////
/*     Fragmentation Model Class     */  
///////////////////////////////////////

#ifndef CHAMBERS

// Ref J.D. Kominami et al. 2018 in preparation
class Collision : public Collision0 {
 public:
    static PS::F64 a_frag;

    PS::S32 collisionOutcome();

    static void readParameter(std::string name,
                              std::string value){
        if ( name == "a_frag" ){
            a_frag = std::atof(value.c_str());
        }
    }
    static void showParameter() {
        std::cout << std::fixed<<std::setprecision(5)
                  << "a_frag        = " << a_frag << std::endl;
                  
    }
    static void showParameter(std::ofstream & fout) {
        fout << std::fixed<<std::setprecision(5)
             << "a_frag        = " << a_frag << std::endl;
             
    }
};

PS::F64 Collision::a_frag = 0.5;

inline PS::S32 Collision::collisionOutcome()
{
    //Number of Fragments
    n_frag = (PS::S32)floor(a_frag * mass_imp / m_min);
    if ( n_frag < 2 ) n_frag = 0;
    if ( n_frag > 10 ) n_frag = 10;
    
    //Mass
    mass_frag = n_frag ? (a_frag * mass_imp / n_frag) : 0.;
    
    //Pos & Vel
    pos_imp_new = pos_tar_new = (mass_imp*pos_imp + mass_tar*pos_tar)/(mass_imp+mass_tar);
    vel_imp_new = vel_tar_new = (mass_imp*vel_imp + mass_tar*vel_tar)/(mass_imp+mass_tar);

    HitAndRun = false;
    
    return n_frag;
}

#else

// Ref J.E. Chambers 2013; S.T. Stewart & Z.M. Leinhardt 2012 
class Collision : public Collision0 {
 public:
    static PS::F64 dens;
    static PS::F64 c_s;
    static PS::F64 mu_;
    static PS::F64 eta_;

    static PS::F64 eps_n;
    static PS::F64 eps_t;

    PS::S32 collisionOutcome();

    static void readParameter(std::string name,
                              std::string value){
        const PS::F64 L_MKS = 149597870700;
        const PS::F64 L_CGS = 14959787070000;
        const PS::F64 M_MKS = 1.9884e30;
        const PS::F64 M_CGS = 1.9884e33;
        if ( name == "dens_imp" ){
            dens = getvalue(value, M_MKS/(L_MKS*L_MKS*L_MKS), M_CGS/(L_CGS*L_CGS*L_CGS));
        } else if ( name == "c_s" ){
            c_s = getvalue(value, 1., 1.);
        } else if ( name == "mu_" ){
            mu_ = getvalue(value, 1., 1.);
        } else if ( name == "eta_" ){
            eta_ = getvalue(value, 1., 1.);
        } else if ( name == "eps_n" ){
            eps_n = getvalue(value, 1., 1.);
        } else if ( name == "eps_t" ){
            eps_t = getvalue(value, 1., 1.);
        } 
    }
    static void showParameter() {
        const PS::F64 L = 149597870700;
        const PS::F64 M = 1.9884e30;
        std::cout << std::scientific << std::setprecision(15)
                  << "dens_imp      = " << dens << "\t(" << dens*M/(L*L*L) << " kg/m^3)"<< std::endl
                  << "c_s           = " << c_s << std::endl
                  << "mu_           = " << mu_ << std::endl
                  << "eta_          = " << eta_ << std::endl
                  << "eps_n         = " << eps_n << std::endl
                  << "eps_t         = " << eps_t << std::endl;
    }
    static void showParameter(std::ofstream & fout) {
        const PS::F64 L = 149597870700;
        const PS::F64 M = 1.9884e30;
        fout << std::scientific << std::setprecision(15)
             << "dens_imp      = " << dens << "\t(" << dens*M/(L*L*L) << " kg/m^3)"<< std::endl
             << "c_s           = " << c_s << std::endl
             << "mu_           = " << mu_ << std::endl
             << "eta_          = " << eta_ << std::endl
             << "eps_n         = " << eps_n << std::endl
             << "eps_t         = " << eps_t << std::endl;
    }
};

PS::F64 Collision::dens  = 1.68373e6;
PS::F64 Collision::c_s   = 1.8;
PS::F64 Collision::mu_   = 1./3.;
PS::F64 Collision::eta_  = -3./2.;
PS::F64 Collision::eps_n = 1.;
PS::F64 Collision::eps_t = 1.;

inline PS::S32 Collision::collisionOutcome()
{
    //const PS::F64 dens = FPGrav::dens;
    const PS::F64 c1 = 2.43;
    const PS::F64 c2 = -0.0408;
    const PS::F64 c3 = 1.86;
    const PS::F64 c4 = 1.08;
    const PS::F64vec ximp = pos_imp - pos_tar;
    const PS::F64vec vimp = vel_imp - vel_tar;
    const PS::F64vec vel_g = (mass_imp*vel_imp + mass_tar*vel_tar)/(mass_imp+mass_tar);

    PS::F64 R_tar = pow(0.75*mass_tar/(M_PI*FPGrav::dens), 1./3.);
    PS::F64 R_imp = pow(0.75*mass_imp/(M_PI*FPGrav::dens), 1./3.);
    PS::F64 R = pow(0.75*(mass_imp + mass_tar)/(M_PI*FPGrav::dens), 1./3.);
    PS::F64 b0 = sin(col_angle);
    PS::F64 b  = b0;
    PS::F64 vel_impact = sqrt(vimp*vimp); // impact velocity corrected for enhancement factor
    PS::F64 vel_escf = sqrt(2.*(mass_imp + mass_tar)/sqrt(ximp*ximp));
    
    // Correction for enhancement factor
#if 1
    if ( f > 1. ) {
        PS::F64 vel_impactf = vel_impact;
        PS::F64 theta2 = vel_escf*vel_escf / (vel_impactf*vel_impactf-vel_escf*vel_escf);
        PS::F64 h = f * sqrt((1.+theta2) / (1.+f*theta2));
      
        vel_impact = sqrt(vel_impactf*vel_impactf + (f-1.)*vel_escf*vel_escf);
        b = f/h * vel_impactf/vel_impact * b;
        if ( b > 1. ) b = 1.;
    }
#endif
    
    PS::F64 gamma = mass_imp / mass_tar;
    PS::F64 Gamma = (1.-gamma)*(1.-gamma)/((1.+gamma)*(1.+gamma));

    PS::F64 mu = mass_imp * mass_tar / ( mass_imp + mass_tar );
    PS::F64 Q = mu * vel_impact*vel_impact / (2.*(mass_imp+mass_tar));
    PS::F64 Q_0 = c_s * 0.8 * M_PI * dens * R*R;

    PS::F64 l = (R_tar + R_imp)*(1.-b);
    PS::F64 alpha = ( l < 2.*R_imp ) ? (3.*R_imp*l*l-l*l*l)/(4.*R_imp*R_imp*R_imp) : 1.;
    PRC(b);PRC(l/R_imp);PRL(alpha);
    //assert ( 0. <= alpha && alpha <= 1. );
    PS::F64 mu_alpha = alpha * mass_imp * mass_tar / ( alpha * mass_imp + mass_tar );
    
    PS::F64 Q_s = pow(mu/mu_alpha, 2.-1.5*mu_) * pow(0.25*(gamma+1.)*(gamma+1.)/gamma, 2./(3.*mu_)-1.) * Q_0;

    PS::F64 mass_rem;
    if ( Q < 1.8 * Q_s ) {
        mass_rem = (mass_imp + mass_tar) * (1. - Q/(2.*Q_s));
    } else {
        mass_rem = 0.1 * (mass_imp + mass_tar) * pow(Q/(1.8*Q_s), eta_);
    }
    if ( mass_rem < m_min ) mass_rem = m_min;
    
    PS::F64 mass_f;
    if ( b <= R_tar/(R_tar + R_imp) || mass_rem < mass_tar ) {
        mass_f = mass_imp + mass_tar - mass_rem;
        pos_imp_new = pos_tar_new = (mass_imp*pos_imp + mass_tar*pos_tar)/(mass_imp+mass_tar);
        vel_imp_new = vel_tar_new = (mass_imp*vel_imp + mass_tar*vel_tar)/(mass_imp+mass_tar);
        HitAndRun = false;
        
        //Number of Fragments
        assert ( mass_f >= 0. );
        n_frag = (PS::S32)floor(mass_f / m_min);
        if ( n_frag < 2 ) n_frag = 0;
        if ( n_frag > 10 ) n_frag = 10;
        mass_frag = n_frag ? (mass_f / n_frag) : 0.;
        
    } else if ( sqrt(vimp*vimp)/vel_escf < (c1*Gamma+c3)*pow(1.-b0, 5./2.) + c2*Gamma+c4 ){
        mass_f = 0.;
        pos_imp_new = pos_tar_new = (mass_imp*pos_imp + mass_tar*pos_tar)/(mass_imp+mass_tar);
        vel_imp_new = vel_tar_new = (mass_imp*vel_imp + mass_tar*vel_tar)/(mass_imp+mass_tar);
        HitAndRun = false;

        //Number of Fragments
        n_frag = 0;
        mass_frag = 0.;
        
    } else {
        PS::F64 beta = ( l < 2.*R_tar ) ? (3.*R_tar*l*l-l*l*l)/(4.*R_tar*R_tar*R_tar) : 1.;
        PRC(l/R_imp);PRL(beta);
        //assert ( 0. <= beta && beta <= 1. );
        gamma = beta * mass_tar / mass_imp;
        R = pow(0.75*(mass_imp + beta*mass_tar)/(M_PI*FPGrav::dens), 1./3.);
        Q_0 = c_s * 0.8 * M_PI * dens * R*R;
        Q_s = pow(0.25*(gamma+1.)*(gamma+1.)/gamma, 2./(3.*mu_)-1.) * Q_0;

        mu = mass_imp * beta*mass_tar / ( mass_imp + beta*mass_tar );
        Q = mu * vel_impact*vel_impact / (2.*(mass_imp + beta*mass_tar));

        if ( Q < 1.8 * Q_s ) {
            mass_rem = (mass_imp + beta*mass_tar) * (1. - Q/(2.*Q_s));
        } else {
            mass_rem = 0.1 * (mass_imp + beta*mass_tar) * pow(Q/(1.8*Q_s), eta_);
        }

        if ( mass_rem < m_min ) mass_rem = m_min;
        mass_f = mass_imp - mass_rem;
        if ( mass_f < 0. ) mass_f = 0.;
        HitAndRun = true;

        //Number of Fragments
        assert ( mass_f >= 0. );
        n_frag = (PS::S32)floor(mass_f / m_min);
        if ( n_frag < 2 ) n_frag = 0;
        if ( n_frag > 10 ) n_frag = 10;
        mass_frag = n_frag ? (mass_f / n_frag) : 0.;

        //Pos & Vel
        pos_imp_new = pos_imp;
        pos_tar_new = pos_tar;
        PS::F64vec e_n = ximp/sqrt(ximp*ximp);
        PS::F64vec e_t = vimp - (vimp*e_n) * e_n;
        e_t /= sqrt(e_t*e_t);
        
        PS::F64 vel_imp_n = (vel_imp-vel_g)*e_n;
        PS::F64 vel_tar_n = (vel_tar-vel_g)*e_n;
        PS::F64 vel_imp_t = (vel_imp-vel_g)*e_t;
        PS::F64 vel_tar_t = (vel_tar-vel_g)*e_t;
        
        PS::F64 vel_imp_new_n = (mass_imp-eps_n*mass_tar)*vel_imp_n + (1.+eps_n)*mass_tar*vel_tar_n;
        PS::F64 vel_tar_new_n = (1.+eps_n)*mass_imp*vel_imp_n + (-eps_n*mass_imp+mass_tar)*vel_tar_n;
        PS::F64 vel_imp_new_t = (mass_imp+eps_t*mass_tar)*vel_imp_t + (1.-eps_t)*mass_tar*vel_tar_t;
        PS::F64 vel_tar_new_t = (1.-eps_t)*mass_imp*vel_imp_t + (eps_t*mass_imp+mass_tar)*vel_tar_t;
        vel_imp_new_n /= (mass_imp+mass_tar);
        vel_tar_new_n /= (mass_imp+mass_tar);
        vel_imp_new_t /= (mass_imp+mass_tar);
        vel_tar_new_t /= (mass_imp+mass_tar);
        
        vel_imp_new = vel_g + mass_imp/(mass_imp - n_frag*mass_frag)*(vel_imp_new_n*e_n + vel_imp_new_t*e_t);
        //vel_imp_new = vel_g + vel_imp_new_n*e_n + vel_imp_new_t*e_t;
        vel_tar_new = vel_g + vel_tar_new_n*e_n + vel_tar_new_t*e_t;
    }
    
    return n_frag;
}
#endif
