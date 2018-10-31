#pragma once

class Collision0{
 public:
    PS::F64 time;

    PS::S32 id_imp;
    PS::S32 id_tar;
    PS::S32 id_c_imp;
    PS::S32 id_c_tar;

    PS::F64vec pos_imp;
    PS::F64vec pos_tar;
    PS::F64vec vel_imp;
    PS::F64vec vel_tar;
    PS::F64vec pos_imp_new;
    PS::F64vec pos_tar_new;
    PS::F64vec vel_imp_new;
    PS::F64vec vel_tar_new;
    PS::F64 col_angle;

    PS::F64 mass_imp;
    PS::F64 mass_tar;
    PS::F64 mass_frag;
    
    PS::F64 edisp;
    PS::F64 edisp_d;
    
    PS::S32 n_frag;
    PS::S32 id_frag;
    PS::S32 id_c_frag;

    //PS::S32 rank;

    bool HitAndRun;

    static PS::F64 f;
    static PS::F64 m_min;

    PS::S32 getNumberOfFragment() const { return n_frag; }
    PS::S32 getFragmentID() const { return id_frag; }
    PS::S32 getFragmentIDCluster() const { return id_c_frag; }
    PS::F64 getEnergyDissipation() const { return edisp; }
    PS::F64 getHardEnergyDissipation() const { return edisp_d; }

    template <class Tpsys>
    void inputPair(Tpsys & pp,
                   std::multimap<PS::S32,PS::S32> & merge_list,
                   std::pair<PS::S32,PS::S32> col_pair);
    PS::S32 collisionOutcome(){
        mass_frag = 0.;

        pos_imp_new = pos_tar_new = (mass_imp*pos_imp + mass_tar*pos_tar)/(mass_imp+mass_tar);
        vel_imp_new = vel_tar_new = (mass_imp*vel_imp + mass_tar*vel_tar)/(mass_imp+mass_tar);
        
        n_frag = id_frag = id_c_frag = 0;
        HitAndRun = false;
        return n_frag;
    }
    template <class Tpsys>
    void setParticle(Tpsys & pp,
                     std::multimap<PS::S32,PS::S32> & merge_list,
                     PS::S32 & id_next);
    
    template <class Tpsys>
    PS::F64 calcEnergyDissipation(Tpsys & pp,
                                  std::multimap<PS::S32,PS::S32> & merge_list);
    template <class Tpsys>
    void setNeighbors(Tpsys & pp);
    
    void setNewFragmentID( std::map<PS::S32, PS::S32> id_old2new ){
        if ( id_imp < 0 ) id_imp = id_old2new.at(id_imp);
        if ( id_tar < 0 ) id_tar = id_old2new.at(id_tar);
        if ( id_frag < 0 ) id_frag = id_old2new.at(id_frag);
    }

    void write2File(std::ofstream & fp) const {
        PS::F64vec ximp = pos_imp - pos_tar;
        PS::F64vec vimp = vel_imp - vel_tar;
        fp << std::fixed<<std::setprecision(8)<< this->time << "\t"
           << this->id_imp << "\t" << this->id_tar  << "\t"
           << this->n_frag << "\t" << this->id_frag << "\t"
           << std::scientific<<std::setprecision(15)
           << this->mass_imp  << "\t" << this->mass_tar  << "\t" << this->mass_frag << "\t"
           << sqrt(ximp*ximp) << "\t" << sqrt(vimp*vimp) << "\t"
           << this->col_angle << "\t" << ((this->HitAndRun) ? 1:0) //<< "\t"
            //<< this->rank
           << std::endl;
    }
};

PS::F64 Collision0::f = 1.;
PS::F64 Collision0::m_min = 9.426627927538057e-12;

template <class Tpsys>
inline void Collision0::inputPair(Tpsys & pp,
                                  std::multimap<PS::S32,PS::S32> & merge_list,
                                  std::pair<PS::S32,PS::S32> col_pair){
    id_c_imp = col_pair.first;
    id_c_tar = col_pair.second;
    using iterator = std::multimap<PS::S32,PS::S32>::iterator;
    std::pair<iterator, iterator> imp_range = merge_list.equal_range(id_c_imp);
    std::pair<iterator, iterator> tar_range = merge_list.equal_range(id_c_tar);
        
    id_imp = pp[id_c_imp].id;
    id_tar = pp[id_c_tar].id;
    time = pp[id_c_imp].time;
    assert ( time == pp[id_c_tar].time );
        
    pos_imp = pp[id_c_imp].pos;
    pos_tar = pp[id_c_tar].pos;
    vel_imp = pp[id_c_imp].vel;
    vel_tar = pp[id_c_tar].vel;
    PS::F64vec ximp = pos_imp - pos_tar;
    PS::F64vec vimp = vel_imp - vel_tar;
    col_angle = acos(-(ximp*vimp)/(sqrt(ximp*ximp)*sqrt(vimp*vimp)));

    mass_imp = pp[id_c_imp].mass;
    mass_tar = pp[id_c_tar].mass;
    if ( pp[id_c_imp].isMerged ) {
        for (iterator it = imp_range.first; it != imp_range.second; ++it){
            PS::S32 id_i = it->second;
            mass_imp += pp[id_i].mass;
            assert( time == pp[id_i].time );
        }
    }
    if ( pp[id_c_tar].isMerged ) {
        for (iterator it = tar_range.first; it != tar_range.second; ++it){
            PS::S32 id_j = it->second;
            mass_tar += pp[id_j].mass;
            assert( time == pp[id_j].time );
        }
    }

    //rank = PS::Comm::getRank();
    PRC(id_imp);PRL(id_tar);
}

template <class Tpsys>
inline void Collision0::setParticle(Tpsys & pp,
                                    std::multimap<PS::S32,PS::S32> & merge_list,
                                    PS::S32 & id_next)
{
    using iterator = std::multimap<PS::S32,PS::S32>::iterator;
    std::pair<iterator, iterator> imp_range = merge_list.equal_range(id_c_imp);
    std::pair<iterator, iterator> tar_range = merge_list.equal_range(id_c_tar);
    const PS::F64 dens  = FPGrav::dens;
    const PS::F64 eps2  = FPGrav::eps2;

    PS::F64vec pos_g = (mass_imp*pos_imp + mass_tar*pos_tar)/(mass_imp+mass_tar);
    PS::F64vec vel_g = (mass_imp*vel_imp + mass_tar*vel_tar)/(mass_imp+mass_tar);
    
    ///////////////////
    /*   Accretion   */
    ///////////////////
    //Mass & Radius
    pp[id_c_imp].mass -= n_frag * mass_frag;
    if ( HitAndRun ) {
        pp[id_c_imp].setRPlanet(mass_imp - n_frag*mass_frag);
    } else {
        pp[id_c_imp].r_planet = 0.;
        pp[id_c_tar].setRPlanet(mass_imp + mass_tar - n_frag*mass_frag);
    }

    // ID
    if ( !HitAndRun ) {
        pp[id_c_imp].id = pp[id_c_tar].id;
        pp[id_c_imp].isDead = true;
        for (iterator it = imp_range.first; it != imp_range.second; ++it){
            PS::S32 id_i = it->second;
            pp[id_i].id = pp[id_c_tar].id;
            pp[id_i].isDead = true;
        }
        pp[id_c_tar].isMerged = true;
    }
    
    //Add Fragments
    PS::F64vec ximp = pos_imp - pos_tar;
    PS::F64vec vimp = vel_imp - vel_tar;
    PS::F64 r_frag = 2. * f * pow(0.75*(mass_imp+mass_tar-n_frag*mass_frag)/(M_PI*dens), 1./3.);
    //if ( HitAndRun ) r_frag = 2. * f * ( pow(0.75*mass_imp/(M_PI*dens), 1./3.) + pow(0.75*mass_tar/(M_PI*dens), 1./3.) );
    PS::F64 r2_frag = r_frag*r_frag + eps2;
    PS::F64 r_frag_inv = sqrt( 1. / r2_frag );
    PS::F64 v_frag = 1.05 * sqrt( 2. * (mass_imp+mass_tar-n_frag*mass_frag) * r_frag_inv );
    //PS::F64vec e0 = vimp;
    //PS::F64vec e1 = ximp - ((ximp*vimp)/(vimp*vimp))*vimp;
    PS::F64vec e0 = ximp;
    PS::F64vec e1 = vimp - ((vimp*ximp)/(ximp*ximp))*ximp;
    e0 = e0 / sqrt(e0*e0);
    e1 = e1 / sqrt(e1*e1);
    
    PS::F64 theta0 = 2. * M_PI * drand48();
    PS::F64vec mv_frag = 0.;
    id_c_frag = n_frag ? pp.size() : 0;
    id_frag   = n_frag ? ( - (100*pp[id_c_imp].id_cluster+id_next+1) ) : 0;
    for ( PS::S32 i=0; i<n_frag; i++ ){
        FPGrav new_frag;
        new_frag.mass = mass_frag;
        new_frag.setRPlanet();

        PS::F64 theta = theta0 + 2.*M_PI*i/n_frag;
        //if ( !HitAndRun ) {
        new_frag.pos = pos_g + r_frag * ( cos(theta)*e0 + sin(theta)*e1 );
        new_frag.vel = vel_g + v_frag * ( cos(theta)*e0 + sin(theta)*e1 );
        //} else {
        //theta = (fmod(theta, 2.*M_PI) - M_PI)/2.;
        //new_frag.pos = pos_imp_new + r_frag * ( cos(theta)*e0 + sin(theta)*e1 );
        //new_frag.vel = vel_imp_new + v_frag * ( cos(theta)*e0 + sin(theta)*e1 );
        //}
        PRC(new_frag.pos);PRC(new_frag.vel);PRC(r_frag);PRL(v_frag);
        mv_frag += new_frag.mass * new_frag.vel;

#ifdef USE_INDIVIDUAL_RADII
        new_frag.r_out    = pp[id_c_imp].r_out;
        new_frag.r_search = pp[id_c_imp].r_search;
#endif
        new_frag.time  = time;
        new_frag.dt    = pp[id_c_imp].dt;
        new_frag.phi   = pp[id_c_imp].phi;
        new_frag.phi_d = pp[id_c_imp].phi_d;
        new_frag.phi_s = pp[id_c_imp].phi_s;

        new_frag.id         = id_frag - i;
        new_frag.id_cluster = pp[id_c_imp].id_cluster;
        new_frag.inDomain = true;
        new_frag.isDead   = false;
        new_frag.isMerged = false;
        id_next ++;
        
        pp.push_back(new_frag);
    }
    if ( n_frag ) assert ( pp.size() == id_c_frag + n_frag );

    //Pos & Vel
    //if ( HitAndRun ) vel_imp_new -= mv_frag / (mass_imp-n_frag*mass_frag);
    pp[id_c_imp].pos = pos_imp_new;
    pp[id_c_imp].vel = vel_imp_new;
    for (iterator it = imp_range.first; it != imp_range.second; ++it){
        PS::S32 id_i = it->second;
        pp[id_i].pos = pos_imp_new;
        pp[id_i].vel = vel_imp_new;
    }    
    pp[id_c_tar].pos = pos_tar_new;
    pp[id_c_tar].vel = vel_tar_new;
    for (iterator it = tar_range.first; it != tar_range.second; ++it){
        PS::S32 id_j = it->second;
        pp[id_j].pos = pos_tar_new;
        pp[id_j].vel = vel_tar_new;
    }
    if ( !HitAndRun ) {
            assert( pos_imp_new == pos_tar_new );
            assert( vel_imp_new == vel_tar_new );
    }
}

template <class Tpsys>
inline PS::F64 Collision0::calcEnergyDissipation(Tpsys & pp,
                                                 std::multimap<PS::S32,PS::S32> & merge_list)
{
    using iterator = std::multimap<PS::S32,PS::S32>::iterator;
    std::pair<iterator, iterator> imp_range = merge_list.equal_range(id_c_imp);
    std::pair<iterator, iterator> tar_range = merge_list.equal_range(id_c_tar);
    PS::F64 mass_i = pp[id_c_imp].mass + n_frag*mass_frag;
    
    const PS::F64 eps2  = EPGrav::eps2;
    const PS::F64 m_sun = FPGrav::m_sun;
    
    ///////////////////////////
    /*   Energy Dissipation  */
    ///////////////////////////
    //Kinetic energy
    PS::F64 e_kin = 0.;
#if 0
    PS::F64vec vel_rel = vel_imp - vel_tar;
    PS::F64vec vel_g = (mass_imp*vel_imp + mass_tar*vel_tar)/(mass_imp+mass_tar);   
    e_kin -= 0.5 * mass_imp*mass_tar/(mass_imp+mass_tar) * vel_rel*vel_rel;

    vel_rel = vel_imp_new - vel_tar_new;
    e_kin += 0.5 * (mass_imp-n_frag*mass_frag)*mass_tar/(mass_imp+mass_tar-n_frag*mass_frag) * vel_rel*vel_rel;
    for ( PS::S32 i=0; i<n_frag; i++ ){
        PS::S32 id_f = id_c_frag + i;
        vel_rel = pp[id_f].vel - vel_g;
        e_kin += 0.5 * pp[id_f].mass * vel_rel*vel_rel;
    }
#else
    e_kin -= mass_imp * vel_imp*vel_imp + mass_tar * vel_tar*vel_tar;

    e_kin += (mass_imp-n_frag*mass_frag) * vel_imp_new*vel_imp_new + mass_tar * vel_tar_new*vel_tar_new;
    for ( PS::S32 i=0; i<n_frag; i++ ){
        PS::S32 id_f = id_c_frag + i;
        e_kin += pp[id_f].mass *pp[id_f].vel*pp[id_f].vel;
    }
    e_kin *= 0.5;
#endif
    
    //Interactive energy
    PS::F64 e_int = 0.;
    PS::F64 e_int_d = 0.;
    PS::F64 dphi = 0;
    
    PS::F64vec dr = pos_imp - pos_tar;
    PS::F64 dr2 = dr*dr + eps2;
    PS::F64 rinv = sqrt(1./dr2);
    PS::F64 rinv_new;
    dphi = mass_i * rinv;
    e_int   += pp[id_c_tar].mass * dphi;
    e_int_d += pp[id_c_tar].mass * dphi
#ifdef USE_INDIVIDUAL_RADII
        * (1.-cutoff_W2(dr2, pp[id_c_imp].r_out, pp[id_c_tar].r_out));
#else
        * (1.-cutoff_W2(dr2, FPGrav::r_CUT));
#endif
    for (iterator it2 = tar_range.first; it2 != tar_range.second; ++it2){
        PS::S32 id_j = it2->second;
        e_int   += pp[id_j].mass * dphi;
        e_int_d += pp[id_j].mass * dphi
#ifdef USE_INDIVIDUAL_RADII
            * (1.-cutoff_W2(dr2, pp[id_c_imp].r_out, pp[id_j].r_out));
#else
            * (1.-cutoff_W2(dr2, FPGrav::r_CUT));
#endif
    }
        
    for (iterator it = imp_range.first; it != imp_range.second; ++it){
        PS::S32 id_i = it->second;
        dphi = pp[id_i].mass * rinv;
        e_int   += pp[id_c_tar].mass * dphi;
        e_int_d += pp[id_c_tar].mass * dphi
#ifdef USE_INDIVIDUAL_RADII
            * (1.-cutoff_W2(dr2, pp[id_i].r_out, pp[id_c_tar].r_out));
#else
            * (1.-cutoff_W2(dr2, FPGrav::r_CUT));
#endif
        for (iterator it2 = tar_range.first; it2 != tar_range.second; ++it2){
            PS::S32 id_j = it2->second;
            e_int   += pp[id_j].mass * dphi;
            e_int_d += pp[id_j].mass * dphi
#ifdef USE_INDIVIDUAL_RADII
                * (1.-cutoff_W2(dr2, pp[id_i].r_out, pp[id_j].r_out));
#else
                * (1.-cutoff_W2(dr2, FPGrav::r_CUT));
#endif
        }
    }

    if ( HitAndRun ){
        dr = pos_imp_new - pos_tar_new;
        dr2 = dr*dr + eps2;
        rinv = sqrt(1./dr2);
        dphi = -pp[id_c_imp].mass * rinv;
        e_int   += pp[id_c_tar].mass * dphi;
        e_int_d += pp[id_c_tar].mass * dphi
#ifdef USE_INDIVIDUAL_RADII
            * (1.-cutoff_W2(dr2, pp[id_c_imp].r_out, pp[id_c_tar].r_out));
#else
            * (1.-cutoff_W2(dr2, FPGrav::r_CUT));
#endif
        for (iterator it2 = tar_range.first; it2 != tar_range.second; ++it2){
            PS::S32 id_j = it2->second;
            e_int   += pp[id_j].mass * dphi;
            e_int_d += pp[id_j].mass * dphi
#ifdef USE_INDIVIDUAL_RADII
                * (1.-cutoff_W2(dr2, pp[id_c_imp].r_out, pp[id_j].r_out));
#else
                * (1.-cutoff_W2(dr2, FPGrav::r_CUT));
#endif
        }
        
        for (iterator it = imp_range.first; it != imp_range.second; ++it){
            PS::S32 id_i = it->second;
            dphi = -pp[id_i].mass * rinv;
            e_int   += pp[id_c_tar].mass * dphi;
            e_int_d += pp[id_c_tar].mass * dphi
#ifdef USE_INDIVIDUAL_RADII
                * (1.-cutoff_W2(dr2, pp[id_i].r_out, pp[id_c_tar].r_out));
#else
                * (1.-cutoff_W2(dr2, FPGrav::r_CUT));
#endif
            for (iterator it2 = tar_range.first; it2 != tar_range.second; ++it2){
                PS::S32 id_j = it2->second;
                e_int   += pp[id_j].mass * dphi;
                e_int_d += pp[id_j].mass * dphi
#ifdef USE_INDIVIDUAL_RADII
                    * (1.-cutoff_W2(dr2, pp[id_i].r_out, pp[id_j].r_out));
#else
                    * (1.-cutoff_W2(dr2, FPGrav::r_CUT));
#endif
            }
        }
    } else {
        assert ( pos_imp_new == pos_tar_new );
    }
    
    for ( PS::S32 i=0; i<n_frag; i++ ){
        PS::S32 id_f0 = id_c_frag + i;
        dr = pp[id_f0].pos - pp[id_c_imp].pos;
        dr2 = dr*dr + eps2;
        rinv_new = sqrt(1./dr2);
        dphi = -pp[id_f0].mass * rinv_new;
        e_int   += pp[id_c_imp].mass * dphi;
        e_int_d += pp[id_c_imp].mass * dphi
#ifdef USE_INDIVIDUAL_RADII
            * (1.-cutoff_W2(dr2, pp[id_c_imp].r_out, pp[id_f0].r_out));
#else
            * (1.-cutoff_W2(dr2, FPGrav::r_CUT));
#endif
        for (iterator it = imp_range.first; it != imp_range.second; ++it){
            PS::S32 id_i = it->second;
            e_int   += pp[id_i].mass * dphi;
            e_int_d += pp[id_i].mass * dphi
#ifdef USE_INDIVIDUAL_RADII
                * (1.-cutoff_W2(dr2, pp[id_i].r_out, pp[id_f0].r_out));
#else
                * (1.-cutoff_W2(dr2, FPGrav::r_CUT));
#endif
        }

        dr = pp[id_f0].pos - pp[id_c_tar].pos;
        dr2 = dr*dr + eps2;
        rinv_new = sqrt(1./dr2);
        dphi = -pp[id_f0].mass * rinv_new;
        e_int   += pp[id_c_tar].mass * dphi;
        e_int_d += pp[id_c_tar].mass * dphi
#ifdef USE_INDIVIDUAL_RADII
            * (1.-cutoff_W2(dr2, pp[id_c_tar].r_out, pp[id_f0].r_out));
#else
            * (1.-cutoff_W2(dr2, FPGrav::r_CUT));
#endif
        for (iterator it = tar_range.first; it != tar_range.second; ++it){
            PS::S32 id_j = it->second;
            e_int   += pp[id_j].mass * dphi;
            e_int_d += pp[id_j].mass * dphi
#ifdef USE_INDIVIDUAL_RADII
                * (1.-cutoff_W2(dr2, pp[id_j].r_out, pp[id_f0].r_out));
#else
                * (1.-cutoff_W2(dr2, FPGrav::r_CUT));
#endif
        }      
        for ( PS::S32 j=0; j<i; j++ ){
            PS::S32 id_f1 = id_c_frag + j;
            dr = pp[id_f0].pos - pp[id_f1].pos;
            dr2 = dr*dr + eps2;
            rinv_new = sqrt(1./dr2);
            dphi = -pp[id_f0].mass * pp[id_f1].mass * rinv_new;
            e_int   += dphi;
            e_int_d += dphi
#ifdef USE_INDIVIDUAL_RADII
                * (1.-cutoff_W2(dr2, pp[id_f0].r_out, pp[id_f1].r_out));
#else
                * (1.-cutoff_W2(dr2, FPGrav::r_CUT));
#endif
        }
    }
    
#if 1
    for ( PS::S32 i=0; i<pp[id_c_imp].neighbor; i++ ){
        PS::S32 id_nei = pp[id_c_imp].n_hard_list.at(i);
        if ( pp[id_nei].id == pp[id_c_imp].id || pp[id_nei].id == pp[id_c_tar].id ) continue;
        PS::F64 m_nei = pp[id_nei].mass;
        
        dr = pos_imp - pp[id_nei].pos;
        dr2 = dr*dr + eps2;
        rinv = sqrt(1./dr2);
        dr = pp[id_c_imp].pos - pp[id_nei].pos;
        dr2 = dr*dr + eps2;
        rinv_new = sqrt(1./dr2);
        dphi = m_nei * ( rinv - rinv_new );
        e_int   += m_nei * ( mass_i * rinv - pp[id_c_imp].mass * rinv_new );
        e_int_d += m_nei * ( mass_i * rinv - pp[id_c_imp].mass * rinv_new )
#ifdef USE_INDIVIDUAL_RADII
            * (1.-cutoff_W2(dr2, pp[id_c_imp].r_out, pp[id_nei].r_out));
#else
            * (1.-cutoff_W2(dr2, FPGrav::r_CUT));
#endif
        for (iterator it = imp_range.first; it != imp_range.second; ++it){
            PS::S32 id_i = it->second;
            e_int   += pp[id_i].mass * dphi;
            e_int_d += pp[id_i].mass * dphi
#ifdef USE_INDIVIDUAL_RADII
                * (1.-cutoff_W2(dr2, pp[id_i].r_out, pp[id_nei].r_out));
#else
                * (1.-cutoff_W2(dr2, FPGrav::r_CUT));
#endif
        }
        for ( PS::S32 j=0; j<n_frag; j++ ){
            PS::S32 id_f = id_c_frag + j;
            dr = pp[id_f].pos - pp[id_nei].pos;
            dr2 = dr*dr + eps2;
            rinv_new = sqrt(1./dr2);
            dphi = -m_nei * pp[id_f].mass * rinv_new;
            e_int   += dphi;
            e_int_d += dphi
#ifdef USE_INDIVIDUAL_RADII
                * (1.-cutoff_W2(dr2, pp[id_f].r_out, pp[id_nei].r_out));
#else
                * (1.-cutoff_W2(dr2, FPGrav::r_CUT));
#endif
        }
    }
    for ( PS::S32 i=0; i<pp[id_c_tar].neighbor; i++ ){
        PS::S32 id_nei = pp[id_c_tar].n_hard_list.at(i);
        if ( pp[id_nei].id == pp[id_c_imp].id || pp[id_nei].id == pp[id_c_tar].id ) continue;
        PS::F64 m_nei = pp[id_nei].mass;
        
        dr = pos_tar - pp[id_nei].pos;
        dr2 = dr*dr + eps2;
        rinv = sqrt(1./dr2);
        dr = pp[id_c_tar].pos - pp[id_nei].pos;
        dr2 = dr*dr + eps2;
        rinv_new = sqrt(1./dr2);
        dphi = m_nei * ( rinv - rinv_new );
        e_int   += pp[id_c_tar].mass * dphi;
        e_int_d += pp[id_c_tar].mass * dphi
#ifdef USE_INDIVIDUAL_RADII
            * (1.-cutoff_W2(dr2, pp[id_c_tar].r_out, pp[id_nei].r_out));
#else
            * (1.-cutoff_W2(dr2, FPGrav::r_CUT));
#endif
        for (iterator it = tar_range.first; it != tar_range.second; ++it){
            PS::S32 id_j = it->second;
            e_int   += pp[id_j].mass * dphi;
            e_int_d += pp[id_j].mass * dphi
#ifdef USE_INDIVIDUAL_RADII
                * (1.-cutoff_W2(dr2, pp[id_j].r_out, pp[id_nei].r_out));
#else
                * (1.-cutoff_W2(dr2, FPGrav::r_CUT));
#endif
        }
    }
#endif
    
    // Gravitational energy
    PS::F64 e_sun = 0.;
    e_sun +=  mass_imp/sqrt(pos_imp*pos_imp+eps2)
        + mass_tar/sqrt(pos_tar*pos_tar+eps2);

    dr = pp[id_c_tar].pos;
    dr2 = dr*dr + eps2;
    rinv_new = sqrt(1./dr2);
    e_sun -= mass_tar * rinv_new;
    
    dr = pp[id_c_imp].pos;
    dr2 = dr*dr + eps2;
    rinv_new = sqrt(1./dr2);
    e_sun -= (mass_imp-n_frag*mass_frag) * rinv_new;

    for ( PS::S32 i=0; i<n_frag; i++ ){
        PS::S32 id_f = id_c_frag + i;
        dr = pp[id_f].pos;
        dr2 = dr*dr + eps2;
        rinv_new = sqrt(1./dr2);
        e_sun -= pp[id_f].mass * rinv_new;
    }
    e_sun *= m_sun;
    
    edisp = e_kin + e_int + e_sun;
    edisp_d = e_kin + e_int_d + e_sun;
    PRC(e_kin); PRC(e_int);PRC(e_int_d); PRC(e_sun); PRL(edisp);

    return edisp;
}

template <class Tpsys>
inline void Collision0::setNeighbors(Tpsys & pp)
{   
    ///////////////////////
    /*   Set Neighbors   */
    ///////////////////////
    for ( PS::S32 i=0; i<n_frag; i++ ){
        PS::S32 id_f = id_c_frag + i;
        pp[id_f].neighbor = pp[id_c_imp].neighbor + n_frag;
        pp[id_f].n_hard_list.clear();
        pp[id_f].n_hard_list.push_back(id_c_imp);
        for ( PS::S32 j=0; j<pp[id_c_imp].neighbor; j++ ){
            PS::S32 id_nei = pp[id_c_imp].n_hard_list.at(j);
            pp[id_f].n_hard_list.push_back(id_nei);
            pp[id_nei].n_hard_list.push_back(id_f);
            pp[id_nei].neighbor ++;
        }
        for ( PS::S32 j=0; j<n_frag; j++ ){
            if ( i != j ) pp[id_f].n_hard_list.push_back(id_c_frag + j);
        }
        assert ( pp[id_f].neighbor == pp[id_f].n_hard_list.size() );
    }

    pp[id_c_imp].neighbor += n_frag;
    for ( PS::S32 i=0; i<n_frag; i++ ){
        pp[id_c_imp].n_hard_list.push_back(id_c_frag + i);
    }
    assert ( pp[id_c_imp].neighbor == pp[id_c_imp].n_hard_list.size() );
    for ( PS::S32 i=0; i<pp[id_c_imp].neighbor; i++ ){
        PS::S32 id_nei = pp[id_c_imp].n_hard_list.at(i);
        assert ( pp[id_nei].neighbor == pp[id_nei].n_hard_list.size() );
    }
}
