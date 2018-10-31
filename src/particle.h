#pragma once

class ForceGrav{
public:
    PS::F64vec acc;
    PS::F64    phi;
    
    void clear(){
        acc = 0.;
        phi = 0.;
    }
};


class EPGrav{
  public:
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64vec acc;

    PS::F64 mass;
    
    PS::S32 id;
    PS::S32 id_local;
    PS::S32 myrank;

#ifdef USE_INDIVIDUAL_RADII
    PS::F64 r_out;
    PS::F64 r_search;
#endif

    static PS::F64 eps2;
    static PS::F64 r_CUT;
    static PS::F64 r_SEARCH;
    static PS::F64 R_cut;
    static PS::F64 R_search0;
    static PS::F64 R_search1;
    static PS::F64 gamma;

    PS::F64 getRSearch() const{
        return r_SEARCH;
    }
    PS::F64vec getPos() const {
        return pos;
    }
    PS::F64 getCharge() const {
        return mass;
    }

    void copyFromFP(const EPGrav & fp){
        pos = fp.pos;
        vel = fp.vel;
        acc = fp.acc;
        
        mass = fp.mass;
        
        id = fp.id;
        id_local = fp.id_local;
        myrank = fp.myrank;
        
#ifdef USE_INDIVIDUAL_RADII
        r_out    = fp.r_out;
        r_search = fp.r_search;
#endif
    }
};

PS::F64 EPGrav::eps2   = 0.;
PS::F64 EPGrav::r_CUT;
PS::F64 EPGrav::r_SEARCH;
PS::F64 EPGrav::R_cut     = 1.;
PS::F64 EPGrav::R_search0 = 1.;
PS::F64 EPGrav::R_search1 = 1.;
PS::F64 EPGrav::gamma     = 0.1;


class FPGrav : public EPGrav {
 public:
    PS::F64vec acc_d;
    PS::F64vec acc_gd;
    PS::F64vec jerk;

    PS::F64 phi_s;
    PS::F64 phi_d;
    PS::F64 phi;
    static PS::F64 m_sun;
    static PS::F64 dens;
    
    PS::F64 time;
    PS::F64 dt;
    static PS::F64 dt_tree;
    static PS::F64 dt_min;
    static PS::F64 eta;
    static PS::F64 eta_0;
    static PS::F64 alpha2;
    
    PS::F64 r_planet;
    static PS::F64 rHill_min;

    PS::S32 id_cluster;
    PS::S32 neighbor;
    
    bool inDomain;
    bool isDead;
    bool isMerged;

    PS::F64 getSemimajorAxis() const {
        return 1.0 / (2.0/sqrt(pos*pos) - vel*vel/m_sun);
    }
    PS::F64 getEccentricity(PS::F64 & ax) const {
        PS::F64 r = sqrt(pos*pos);
        PS::F64 rv = pos*vel;
        ax = 1.0 / (2.0/r - vel*vel/m_sun);
        PS::F64 ecccosu = 1. - r/ax;
        PS::F64 eccsinu = rv/sqrt(m_sun*ax);
        return sqrt(ecccosu*ecccosu + eccsinu*eccsinu);
    }
    PS::F64 getInclination(PS::F64vec & h) const {
        h.x = pos.y*vel.z - pos.z*vel.y;
        h.y = pos.z*vel.x - pos.x*vel.z;
        h.z = pos.x*vel.y - pos.y*vel.x;
        return atan2(sqrt(h.x*h.x + h.y*h.y), h.z);
    }
    PS::F64 getRHill() const {
        PS::F64 ax = 1.0 / (2.0/sqrt(pos*pos) - vel*vel/m_sun);
        return pow(mass/(3.*m_sun), 1./3.) * ax;
    }
#ifdef USE_INDIVIDUAL_RADII
    PS::F64 setMyROutRSearch(PS::F64 vdisp){
        PS::F64 ax = 1.0 / (2.0/sqrt(pos*pos) - vel*vel/m_sun);
        PS::F64 rHill = std::max(pow(mass/(3.*m_sun), 1./3.) * ax, rHill_min);
        
        r_out    = R_cut    *rHill;
        r_search = R_search0*rHill + R_search1*vdisp*dt_tree;
        
        return rHill;
    }
#endif
    static void setRCutRSearch(PS::F64 rHill_glb,
                               PS::F64 vdisp){
        r_CUT    = R_cut    *rHill_glb;
        r_SEARCH = R_search0*rHill_glb + R_search1*vdisp*dt_tree;
    }
    void setRPlanet(PS::F64 m) {
        r_planet = pow(0.75*m/(M_PI*dens), 1./3.);
    }
    void setRPlanet() {
        r_planet = pow(0.75*mass/(M_PI*dens), 1./3.);
    }
    
    void copyFromForce(const ForceGrav & force){
        acc = force.acc;
        phi = force.phi;
    }

    void writeAscii(FILE* fp) const {
        fprintf(fp, "%d\t%20.15e\t%20.15e\t%20.15e\t%20.15e\t%20.15e\t%20.15e\t%20.15e\t%d\n",
        //fprintf(fp, "%d\t%20.15e\t%20.15e\t%20.15e\t%20.15e\t%20.15e\t%20.15e\t%20.15e\t%20.15e\t%20.15e\n", 
                this->id, this->mass,
                this->pos.x, this->pos.y, this->pos.z,
                this->vel.x, this->vel.y, this->vel.z,
                this->neighbor);
        //this->r_out, this->r_search);
    }
    void readAscii(FILE* fp) {
        fscanf(fp, "%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%d\n",
        //fscanf(fp, "%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
               &this->id, &this->mass,
               &this->pos.x, &this->pos.y, &this->pos.z,
               &this->vel.x, &this->vel.y, &this->vel.z,
               &this->neighbor);
               //&this->r_out, &this->r_search);
    }

    void velKick(){
        vel += 0.5*dt_tree*acc;
    }
    //void velKickGas(){
    //    vel += 0.5*dt_tree*acc_d;
    //}
    
    void calcDeltatInitial(){
        PS::F64 dt_next = 0.5*dt_tree;
        
#if 0
        PS::F64 a00 = alpha2 * mass / (r_out * r_out);
#else
        PS::F64 a00 = alpha2 * mass / (r_CUT * r_CUT);
#endif
        PS::F64 dt_2 = (acc_d*acc_d + a00*a00)/(jerk*jerk);
        PS::F64 dt_1 = eta_0 * sqrt(dt_2);
        PS::F64 rem = fmod(time, dt_next);
    
        while( rem != 0.0 ){
            dt_next *= 0.5;
            rem = fmod(time, dt_next);
        }
        while( dt_1 < dt_next ) dt_next *= 0.5;
        if( dt_next < 2.*dt_min ) dt_next = dt_min;
        
        dt = dt_next;
    }
};

PS::F64 FPGrav::m_sun     = 1.;
PS::F64 FPGrav::dens      = 5.049667e6;
PS::F64 FPGrav::dt_tree   = pow2(-5);
PS::F64 FPGrav::dt_min    = pow2(-13);
PS::F64 FPGrav::eta       = 0.01;
PS::F64 FPGrav::eta_0     = 0.001;
PS::F64 FPGrav::alpha2    = 0.01;
PS::F64 FPGrav::rHill_min = 3.15557400894e-04;


class FPHard : public FPGrav {
 public:
    PS::F64vec x0;
    PS::F64vec v0;
    PS::F64vec a0;
    PS::F64vec j0;
    PS::F64vec xp;
    PS::F64vec vp;
    PS::F64vec a2;
    PS::F64vec a3;

    std::vector<PS::S32> n_list;
    std::vector<PS::S32> n_hard_list;
    
    void clearList(){
        std::vector<PS::S32> tmp0, tmp1;
        tmp0.swap(n_list);
        tmp1.swap(n_hard_list);
    }
    void copyList(std::vector<PS::S32> list){
        n_list.resize(list.size());
        std::copy(list.begin(),list.end(),n_list.begin());
    }
    void copyList(PS::S32 * list){
        n_list.clear();
        n_list.reserve(neighbor);
        for ( PS::S32 i=0; i<neighbor; i++ ) n_list.push_back(list[i]);
    }
    void copyHardList(std::vector<PS::S32> list){
        n_hard_list.resize(list.size());
        std::copy(list.begin(),list.end(),n_hard_list.begin());
    }
    void copyHardList(PS::S32 * list){
        n_hard_list.clear();
        n_hard_list.reserve(neighbor);
        for ( PS::S32 i=0; i<neighbor; i++ ) n_hard_list.push_back(list[i]);
    }
    FPHard(){
        jerk = 0.;
        x0 = v0 = a0 = j0 = 0.;
        xp = vp = 0.;
        a2 = a3 = 0.;
        
        clearList();
    }
    FPHard(const FPHard & fp) : FPGrav(fp){
        jerk = fp.jerk;
        x0 = fp.x0;
        v0 = fp.v0;
        a0 = fp.a0;
        j0 = fp.j0;
        xp = fp.xp;
        vp = fp.vp;
        a2 = fp.a2;
        a3 = fp.a3;

        copyList(fp.n_list);
        copyHardList(fp.n_hard_list);
    }
    FPHard(const FPGrav & fp) : FPGrav(fp){
        jerk = 0.;
        x0 = v0 = a0 = j0 = 0.;
        xp = vp = 0.;
        a2 = a3 = 0.;
        
        clearList();
    }
    FPHard &operator=(const FPHard & fp){
        FPGrav::operator=(fp);
        if ( this != &fp ){
            jerk = fp.jerk;
            x0 = fp.x0;
            v0 = fp.v0;
            a0 = fp.a0;
            j0 = fp.j0;
            xp = fp.xp;
            vp = fp.vp;
            a2 = fp.a2;
            a3 = fp.a3;
            
            copyList(fp.n_list);
            copyHardList(fp.n_hard_list);
        }
        return *this;
    }
    FPHard &operator=(const FPGrav & fp){
        FPGrav::operator=(fp);
        if ( this != &fp ){
            jerk = 0.;
            x0 = v0 = a0 = j0 = 0.;
            xp = vp = 0.;
            a2 = a3 = 0.;

            clearList();
        }
        return *this;
    }

    void predict(PS::F64 Dt){
        x0 = pos;
        v0 = vel;
        a0 = acc_d;
        j0 = jerk;
        xp = pos + Dt*vel + Dt*Dt*0.5*acc_d + Dt*Dt*Dt*jerk/6.0;
        vp = vel + Dt*acc_d + Dt*Dt*0.5*jerk;
    }
    void correct(PS::F64 Dt){
        PS::F64 Dt2 = Dt*Dt;
        PS::F64 Dt3 = Dt*Dt*Dt;
        PS::F64 Dt2inv = 1.0/Dt2;
        PS::F64 Dt3inv = 1.0/Dt3;
        a2 = -Dt2inv * (6.0*(a0-acc_d) + Dt*(4.0*j0+2.0*jerk));
        a3 = Dt3inv * (12.0*(a0-acc_d) + 6.0*Dt*(j0+jerk));
        pos = xp + a2*Dt2*Dt2/24.0 + 7.0*a3*Dt2*Dt3/720.0;
        vel = vp + a2*Dt3/6.0 + a3*Dt2*Dt2/24.0;
    }

    void calcDeltatInitial(){
        PS::F64 dt_next = 0.5*dt_tree;

#if 0
        PS::F64 a00 = alpha2 * mass / (r_out * r_out);
#else
        PS::F64 a00 = alpha2 * mass / (r_CUT * r_CUT);
#endif
        PS::F64 dt_2 = (acc_d*acc_d + a00*a00)/(jerk*jerk);
        PS::F64 dt_1 = eta_0 * sqrt(dt_2);
        PS::F64 rem = fmod(time, dt_next);
    
        while( rem != 0.0 ){
            dt_next *= 0.5;
            rem = fmod(time, dt_next);
        }
        while( 2.*dt < dt_next ) dt_next *= 0.5;
        while( dt_1 < dt_next ) dt_next *= 0.5;
        if( dt_next < 2.*dt_min ) dt_next = dt_min;
        
        dt = dt_next;
    }
    void calcDeltat(){
        PS::F64 dt_next = 2.*dt;

#if 0
        PS::F64 a00 = alpha2 * mass / (r_out * r_out);
#else
        PS::F64 a00 = alpha2 * mass / (r_CUT * r_CUT);
#endif
        PS::F64vec a1_2 = a2 + a3*dt;
        PS::F64vec a1_3 = a3;
        PS::F64 b1 = sqrt(acc_d*acc_d + a00*a00);
        PS::F64 b1_dot2 = jerk*jerk;
        PS::F64 b1_dot  = sqrt(b1_dot2);
        PS::F64 b1_22 = a1_2*a1_2;
        PS::F64 b1_2  = sqrt(b1_22);
        PS::F64 b1_3 = sqrt(a1_3*a1_3);
        PS::F64 dt_2 = (b1*b1_2 + b1_dot2)/(b1_dot*b1_3 + b1_22);
        PS::F64 dt_1 = eta * sqrt(dt_2);
        PS::F64 rem = fmod(time, dt_next);
        while(rem != 0.0){
            dt_next *= 0.5;
            rem = fmod(time, dt_next);
        }
        while(dt_1 < dt_next && 0.5*dt < dt_next) dt_next *= 0.5;
        if( dt_next < 2.*dt_min ) dt_next = dt_min;
        
        dt = dt_next;
    }
    
};


template <class Tpsys>
void calcVelDisp(Tpsys & pp,
                   const PS::S32 n_tot,
                   const PS::S32 n_loc,
                   PS::F64 & v_disp)
{
    PS::F64 m_sun   = FPGrav::m_sun;
    PS::F64 ecc_rms_loc = 0.0;
    PS::F64 inc_rms_loc = 0.0;
    PS::F64 v_kep_max_loc = 0.0;
    v_disp = 0.0;

#pragma omp parallel for
    for(PS::S32 i=0; i<n_loc; i++){
        PS::F64 ax;
        PS::F64 ecc = pp[i].getEccentricity(ax);
        PS::F64vec h;
        PS::F64 inc = pp[i].getInclination(h);
        PS::F64 v_kep = sqrt(m_sun/ax);

#pragma omp critical
        {
            ecc_rms_loc += ecc*ecc;
            inc_rms_loc += inc*inc;
            if ( v_kep > v_kep_max_loc ) v_kep_max_loc = v_kep;
        }
    }
    
    PS::F64 ecc_rms2 = PS::Comm::getSum(ecc_rms_loc)/n_tot;
    PS::F64 inc_rms2 = PS::Comm::getSum(inc_rms_loc)/n_tot;
    PS::F64 v_kep_max = PS::Comm::getMaxValue(v_kep_max_loc);
    
    v_disp = (ecc_rms2 + inc_rms2)*v_kep_max;
}

template <class Tpsys>
void setCutoffRadii(Tpsys & pp)
{
    const PS::S32 n_loc = pp.getNumberOfParticleLocal();
    const PS::S32 n_tot = pp.getNumberOfParticleGlobal();
    PS::F64 v_disp;
    //PS::F64 dens = FPGrav::dens;
    
    calcVelDisp(pp, n_tot, n_loc, v_disp);

    PS::F64 rHill_loc = 0.;
#pragma omp parallel for
    for(PS::S32 i=0;i<n_loc;i++){
#ifdef USE_INDIVIDUAL_RADII
        PS::F64 rHill = pp[i].setMyROutRSearch(v_disp);
#else
        PS::F64 rHill = std::max(pp[i].getRHill(), FPGrav::rHill_min);
#endif
        //pp[i].r_planet = pow(0.75*pp[i].mass/(M_PI*dens), 1./3.);
        pp[i].setRPlanet();
#pragma omp critical
        if ( rHill > rHill_loc ) rHill_loc = rHill;
    }
    PS::F64 rHill_glb = PS::Comm::getMaxValue(rHill_loc);
    FPGrav::setRCutRSearch(rHill_glb, v_disp);
    
    //if ( PS::Comm::getRank() == 0) { PRC(EPGrav::r_CUT);PRL(EPGrav::r_SEARCH); }
}
