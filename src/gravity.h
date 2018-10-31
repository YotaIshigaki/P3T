#pragma once

void clearNeighborList(std::vector<std::vector<PS::S32> > & n_list)
{
    const PS::S32 n_loc = n_list.size();
    for(PS::S32 i=0; i<n_loc; i++){
        std::vector<PS::S32> tmp;
        tmp.swap(n_list.at(i));
    }
}

template <class Tpsys>
void velKick(Tpsys & pp){
    PS::S32 n = pp.getNumberOfParticleLocal();
#pragma omp parallel for
    for(PS::S32 i=0; i<n; i++){
        pp[i].velKick();
    }
}

template <class Tpsys>
void calcStarGravity(Tpsys & pp)
{
    PS::F64 eps2 = EPGrav::eps2;
    PS::F64 m_sun = FPGrav::m_sun;
    
    PS::F64vec posi = pp.getPos();
    PS::F64vec veli = pp.vel;
    
    PS::F64vec dr  = - posi;
    PS::F64vec dv  = - veli;
    PS::F64    r2inv = 1.0 / (dr * dr + eps2);
    PS::F64    rinv  = sqrt(r2inv);
    PS::F64    r3inv = rinv * r2inv;
    PS::F64    r5inv = r3inv * r2inv;

    pp.phi_s  = - m_sun * rinv;
    pp.acc_d += m_sun * r3inv * dr;
    pp.jerk  += m_sun * (r3inv*dv -3.0*r5inv*(dr*dv)*dr);
}

template <class Tpsys>
void calcStarGravity_p(Tpsys & pp)
{
    PS::F64 eps2 = EPGrav::eps2;
    PS::F64 m_sun = FPGrav::m_sun;

    PS::F64vec posi = pp.xp;
    PS::F64vec veli = pp.vp;

    PS::F64vec dr  = - posi;
    PS::F64vec dv  = - veli;
    PS::F64    r2inv = 1.0 / (dr * dr + eps2);
    PS::F64    rinv  = sqrt(r2inv);
    PS::F64    r3inv = rinv * r2inv;
    PS::F64    r5inv = r3inv * r2inv;

    pp.acc_d += m_sun * r3inv * dr;
    pp.jerk  += m_sun * (r3inv*dv -3.0*r5inv*(dr*dv)*dr);
}

template <class Tpsys>
void calcStarGravity_c(Tpsys & pp)
{
    PS::F64 eps2 = EPGrav::eps2;
    PS::F64 m_sun = FPGrav::m_sun;

    PS::F64vec posi = pp.getPos();
    PS::F64vec veli = pp.vel;

    PS::F64vec dr  = - posi;
    PS::F64vec dv  = - veli;
    PS::F64    r2inv = 1.0 / (dr * dr + eps2);
    PS::F64    rinv  = sqrt(r2inv);
    PS::F64    r3inv = rinv * r2inv;
    PS::F64    r5inv = r3inv * r2inv;

    pp.acc_d += m_sun * r3inv * dr;
    pp.jerk  += m_sun * (r3inv*dv -3.0*r5inv*(dr*dv)*dr);
}

#if 0
template <class Tpsys>
void calcStarAcc(Tpsys & pp)
{
    PS::F64 eps2 = EPGrav::eps2;
    PS::F64 m_sun = FPGrav::m_sun;

    PS::F64vec posi = pp.getPos();

    PS::F64vec dr  = - posi;
    PS::F64    r2inv = 1.0 / (dr * dr + eps2);
    PS::F64    rinv  = sqrt(r2inv);
    PS::F64    r3inv = rinv * r2inv;

    pp.phi_s  = - m_sun * rinv;
    pp.acc_d += m_sun * r3inv * dr;
}
#endif

template <class Tpsys>
void calcStarJerk(Tpsys & pp)
{
    PS::F64 eps2 = EPGrav::eps2;
    PS::F64 m_sun = FPGrav::m_sun;

    PS::F64vec posi = pp.getPos();
    PS::F64vec veli = pp.vel;

    PS::F64vec dr  = - posi;
    PS::F64vec dv  = - veli;
    PS::F64    r2inv = 1.0 / (dr * dr + eps2);
    PS::F64    rinv  = sqrt(r2inv);
    PS::F64    r3inv = rinv * r2inv;
    PS::F64    r5inv = r3inv * r2inv;

    pp.jerk  += m_sun * (r3inv*dv -3.0*r5inv*(dr*dv)*dr);
}

template <class Tp, class Tpsys>
void calcGravity(Tp & pi,
                 Tpsys & pp)
{
    assert( pi.neighbor != 0 );
    
    PS::F64 eps2  = EPGrav::eps2;
    
    pi.phi_d = 0.0;
    pi.phi_s = 0.0;
    pi.acc_d = 0.0;
    pi.jerk  = 0.0;
    
    calcStarGravity(pi);
    
    PS::S32 pj_id = 0;
    PS::F64vec xi = pi.pos;
    PS::F64vec vi = pi.vel;
    PS::F64vec acci = 0.0;
    PS::F64    phii = 0.0;
    PS::F64vec jerki= 0.0;

#ifndef FORDEBUG
    for(PS::S32 j=0; j<pi.neighbor; j++){
        pj_id = pi.n_hard_list.at(j);
#else
    for(PS::S32 j=0; j<pp.size(); j++){
        pj_id = j;
#endif        
        if ( pi.id == pp[pj_id].id ) continue;
        
        PS::F64vec xj = pp[pj_id].pos;
        PS::F64vec dr  = xj - xi;
        PS::F64 dr2 = dr * dr;
        assert( dr2 != 0.0 );
        dr2 += eps2;
        
        PS::F64vec vj   = pp[pj_id].vel;
        PS::F64vec dv    = vj - vi;
        PS::F64    massj = pp[pj_id].getCharge();

        PS::F64 rij   = sqrt(dr2);
        PS::F64 drdv  = dr*dv;
        PS::F64 r2inv = 1.0 / dr2;
        PS::F64 rinv  = 1./rij;
        PS::F64 r3inv = r2inv * rinv;
        PS::F64 r5inv = r3inv * r2inv;

#ifdef USE_INDIVIDUAL_RADII
        PS::F64 r_out = mutualRadius(pi.r_out, pp[pj_id].r_out);
#else
        PS::F64 r_out = FPGrav::r_CUT;
#endif
        PS::F64 W  = cutoff_W(rij, r_out);
        PS::F64 K  = cutoff_K(rij, r_out);
        PS::F64 dK = cutoff_dK(rij, drdv, r_out);

        phii  -= massj * rinv * (1.-W);
        acci  += massj * r3inv * dr * (1.-K);
        jerki += massj * ( (r3inv*dv -3.*r5inv*(drdv)*dr)*(1.-K) - r3inv*dr*dK );
    }
    pi.acc_d += acci;
    pi.phi_d += phii;
    pi.jerk  += jerki;
}

template <class Tp, class Tpsys>
void calcGravity_p(Tp & pi,
                   Tpsys & pp)
{
    assert( pi.neighbor != 0 );

    PS::F64 eps2  = EPGrav::eps2;
  
    pi.acc_d = 0.0;
    pi.jerk  = 0.0;
    calcStarGravity_p(pi);
        
    PS::S32 pj_id = 0;
    PS::F64vec xpi = pi.xp;
    PS::F64vec vpi = pi.vp;
    PS::F64vec acci = 0.0;
    PS::F64vec jerki= 0.0;

#ifndef FORDEBUG
    for(PS::S32 j=0; j<pi.neighbor; j++){
        pj_id = pi.n_hard_list.at(j);
#else
    for(PS::S32 j=0; j<pp.size(); j++){
        pj_id = j;
#endif
        if ( pi.id == pp[pj_id].id ) continue;

        PS::F64vec xpj = pp[pj_id].xp;
        PS::F64vec dr  = xpj - xpi;
        PS::F64 dr2 = dr * dr;
        assert( dr2 != 0.0 );
        dr2 += eps2;

        PS::F64vec vpj   = pp[pj_id].vp;
        PS::F64vec dv    = vpj - vpi;
        PS::F64    massj = pp[pj_id].getCharge();

        PS::F64 rij   = sqrt(dr2);
        PS::F64 drdv  = dr*dv;
        PS::F64 r2inv = 1.0 / dr2;
        PS::F64 rinv  = 1./rij;
        PS::F64 r3inv = r2inv * rinv;
        PS::F64 r5inv = r3inv * r2inv;

#ifdef USE_INDIVIDUAL_RADII
        PS::F64 r_out = mutualRadius(pi.r_out, pp[pj_id].r_out);
#else
        PS::F64 r_out = FPGrav::r_CUT;
#endif
        PS::F64 K  = cutoff_K(rij, r_out);
        PS::F64 dK = cutoff_dK(rij, drdv, r_out);
        
        acci  += massj * r3inv * dr * (1.-K);
        jerki += massj * ( (r3inv*dv -3.*r5inv*(drdv)*dr)*(1.-K) - r3inv*dr*dK );
    }
    pi.acc_d += acci;
    pi.jerk  += jerki;
}

template <class Tp, class Tpsys>
void calcGravity_c(Tp & pi,
                   Tpsys & pp)
{
    assert( pi.neighbor != 0 );
    
    PS::F64 eps2  = EPGrav::eps2;
    pi.acc_d = 0.0;
    pi.jerk  = 0.0;
    
    calcStarGravity_c(pi);
    
    PS::S32 pj_id = 0;
    PS::F64vec xi = pi.pos;
    PS::F64vec vi = pi.vel;
    PS::F64vec acci = 0.0;
    PS::F64vec jerki= 0.0;

#ifndef FORDEBUG
    for(PS::S32 j=0; j<pi.neighbor; j++){
        pj_id = pi.n_hard_list.at(j);
#else
    for(PS::S32 j=0; j<pp.size(); j++){
        pj_id = j;
#endif
        if ( pi.id == pp[pj_id].id ) continue;
        
        PS::F64vec xpj = pp[pj_id].xp;
        PS::F64vec dr  = xpj - xi;
        PS::F64 dr2 = dr * dr;
        assert( dr2 != 0.0 );
        dr2 += eps2;
        
        PS::F64vec vpj   = pp[pj_id].vp;
        PS::F64vec dv    = vpj - vi;
        PS::F64    massj = pp[pj_id].getCharge();

        PS::F64 rij   = sqrt(dr2);
        PS::F64 drdv  = dr*dv;
        PS::F64 r2inv = 1.0 / dr2;
        PS::F64 rinv  = 1./rij;
        PS::F64 r3inv = r2inv * rinv;
        PS::F64 r5inv = r3inv * r2inv;
        
#ifdef USE_INDIVIDUAL_RADII
        PS::F64 r_out = mutualRadius(pi.r_out, pp[pj_id].r_out);
#else
        PS::F64 r_out = FPGrav::r_CUT;
#endif
        PS::F64 K  = cutoff_K(rij, r_out);
        PS::F64 dK = cutoff_dK(rij, drdv, r_out);
        
        acci  += massj * r3inv * dr * (1.-K);
        jerki += massj * ( (r3inv*dv -3.*r5inv*(drdv)*dr)*(1.-K) - r3inv*dr*dK );
    }
    pi.acc_d += acci;
    pi.jerk  += jerki;
}

template <class Tp, class Tpsys>
void calcJerk(Tp & pi,
              Tpsys & pp)
{   
    PS::F64 eps2  = EPGrav::eps2;
    pi.jerk  = 0.0;
    
    calcStarJerk(pi);
    
    PS::S32 pj_id = 0;
    PS::F64vec xi = pi.pos;
    PS::F64vec vi = pi.vel;
    PS::F64vec jerki= 0.0;

#ifndef FORDEBUG
    for(PS::S32 j=0; j<pi.neighbor; j++){
        pj_id = pi.n_hard_list.at(j);
#else
    for(PS::S32 j=0; j<pp.size(); j++){
        pj_id = j;
#endif
        if ( pi.id == pp[pj_id].id ) continue;
        
        PS::F64vec xj = pp[pj_id].pos;
        PS::F64vec dr  = xj - xi;
        PS::F64 dr2 = dr * dr;
        assert( dr2 != 0.0 );
        dr2 += eps2;
        
        PS::F64vec vj   = pp[pj_id].vel;
        PS::F64vec dv    = vj - vi;
        PS::F64    massj = pp[pj_id].getCharge();

        PS::F64 rij   = sqrt(dr2);
        PS::F64 drdv  = dr*dv;
        PS::F64 r2inv = 1.0 / dr2;
        PS::F64 rinv  = 1./rij;
        PS::F64 r3inv = r2inv * rinv;
        PS::F64 r5inv = r3inv * r2inv;

#ifdef USE_INDIVIDUAL_RADII
        PS::F64 r_out = mutualRadius(pi.r_out, pp[pj_id].r_out);
#else
        PS::F64 r_out = FPGrav::r_CUT;
#endif
        PS::F64 K  = cutoff_K(rij, r_out);
        PS::F64 dK = cutoff_dK(rij, drdv, r_out);

        jerki += massj * ( (r3inv*dv -3.*r5inv*(drdv)*dr)*(1.-K) - r3inv*dr*dK );
    }
    pi.jerk  += jerki;
}


template <class TParticleJ>
void CalcForceLongEP(const EPGrav * ep_i,
                     const PS::S32 n_ip,
                     const TParticleJ * ep_j,
                     const PS::S32 n_jp,
                     ForceGrav * force) {

    PS::F64 eps2 = EPGrav::eps2;
    v4df (*rcp)(v4df) = v4df::rcp_4th;
    PS::S32 nvector = v4df::getVectorLength();

    v4df v0(0.);
    v4df v1(1.);
    v4df rin(EPGrav::gamma*EPGrav::r_CUT);
    v4df rout_rin_inv(1./((1.-EPGrav::gamma)*EPGrav::r_CUT));

#pragma omp parallel for 
    for(PS::S32 i = 0; i < n_ip; i += nvector) {
        PS::F64 buf_px[nvector] __attribute__((aligned(32)));
        PS::F64 buf_py[nvector] __attribute__((aligned(32)));
        PS::F64 buf_pz[nvector] __attribute__((aligned(32)));
        PS::F64 buf_e2[nvector] __attribute__((aligned(32)));
        const PS::S32 nii = std::min(n_ip - i, nvector);
        for(PS::S32 ii = 0; ii < nii; ii++) {
            buf_px[ii] = ep_i[i+ii].pos[0];
            buf_py[ii] = ep_i[i+ii].pos[1];
            buf_pz[ii] = ep_i[i+ii].pos[2];
            buf_e2[ii] = eps2;
        }
        v4df px_i;
        v4df py_i;
        v4df pz_i;
        v4df e2_i;
        px_i.load(buf_px);
        py_i.load(buf_py);
        pz_i.load(buf_pz);
        e2_i.load(buf_e2);
        v4df ax_i(0.0);
        v4df ay_i(0.0);
        v4df az_i(0.0);
        v4df pt_i(0.0);
                
        for(PS::S32 j = 0; j < n_jp; j++) {
            v4df dpx_ij = px_i - v4df(ep_j[j].pos[0]);
            v4df dpy_ij = py_i - v4df(ep_j[j].pos[1]);
            v4df dpz_ij = pz_i - v4df(ep_j[j].pos[2]);

            v4df dr2 = dpx_ij * dpx_ij + dpy_ij * dpy_ij + dpz_ij * dpz_ij + e2_i;

            v4df dr = v4df::max(v4df::sqrt(dr2), rin);
            v4df K = v4df::min( (dr - rin)*rout_rin_inv, v1 );
            
            v4df dr_inv_nK = rcp(dr);
            v4df dr_inv    = rcp(dr) * K;
            v4df dr3_inv = dr_inv_nK * dr_inv_nK * dr_inv;
            

            v4df m_j    = v4df(ep_j[j].mass);
            v4df dg2_ij = m_j * dr3_inv;

            pt_i -= m_j    * dr_inv;
            ax_i -= dpx_ij * dg2_ij;
            ay_i -= dpy_ij * dg2_ij;
            az_i -= dpz_ij * dg2_ij;
        }

        PS::F64 buf_ax[nvector] __attribute__((aligned(32)));
        PS::F64 buf_ay[nvector] __attribute__((aligned(32)));
        PS::F64 buf_az[nvector] __attribute__((aligned(32)));
        PS::F64 buf_pt[nvector] __attribute__((aligned(32)));

        ax_i.store(buf_ax);
        ay_i.store(buf_ay);
        az_i.store(buf_az);
        pt_i.store(buf_pt);

        for(PS::S32 ii = 0; ii < nii; ii++) {
            force[i+ii].acc[0] += buf_ax[ii];
            force[i+ii].acc[1] += buf_ay[ii];
            force[i+ii].acc[2] += buf_az[ii];
            force[i+ii].phi    += buf_pt[ii];
        }
    }
}
 
template <class TParticleJ>
void CalcForceLongSP(const EPGrav * ep_i,
                     const PS::S32 n_ip,
                     const TParticleJ * ep_j,
                     const PS::S32 n_jp,
                     ForceGrav * force) {

    PS::F64 eps2 = EPGrav::eps2;
    v4df (*rcp)(v4df) = v4df::rcp_4th;
    PS::S32 nvector = v4df::getVectorLength();

#ifdef USE_QUAD
    v4df v5(5.0);
    v4df mv2(-2.0);
    v4df v0p5(0.5);
    v4df v1p5(1.5);
#endif
    v4df v0(0.);
    v4df v1(1.);
    v4df rin(EPGrav::gamma*EPGrav::r_CUT);
    v4df rout_rin_inv(1./((1.-EPGrav::gamma)*EPGrav::r_CUT));
    
#pragma omp parallel for 
    for(PS::S32 i = 0; i < n_ip; i += nvector) {
        PS::F64 buf_px[nvector] __attribute__((aligned(32)));
        PS::F64 buf_py[nvector] __attribute__((aligned(32)));
        PS::F64 buf_pz[nvector] __attribute__((aligned(32)));
        PS::F64 buf_e2[nvector] __attribute__((aligned(32)));
        const PS::S32 nii = std::min(n_ip - i, nvector);
        for(PS::S32 ii = 0; ii < nii; ii++) {
            buf_px[ii] = ep_i[i+ii].pos[0];
            buf_py[ii] = ep_i[i+ii].pos[1];
            buf_pz[ii] = ep_i[i+ii].pos[2];
            buf_e2[ii] = eps2;
        }
        v4df px_i;
        v4df py_i;
        v4df pz_i;
        v4df e2_i;
        px_i.load(buf_px);
        py_i.load(buf_py);
        pz_i.load(buf_pz);
        e2_i.load(buf_e2);
        v4df ax_i(0.0);
        v4df ay_i(0.0);
        v4df az_i(0.0);
        v4df pt_i(0.0);
                
        for(PS::S32 j = 0; j < n_jp; j++) {
            v4df dpx_ij = px_i - v4df(ep_j[j].pos[0]);
            v4df dpy_ij = py_i - v4df(ep_j[j].pos[1]);
            v4df dpz_ij = pz_i - v4df(ep_j[j].pos[2]);

            v4df dr2 = dpx_ij * dpx_ij + dpy_ij * dpy_ij + dpz_ij * dpz_ij + e2_i;

            v4df dr = v4df::max(v4df::sqrt(dr2), rin);
            v4df K = v4df::min( (dr - rin)*rout_rin_inv, v1 );

            v4df dr_inv_nK = rcp(dr);
            v4df dr_inv    = rcp(dr) * K;
            v4df dr2_inv = dr_inv_nK * dr_inv_nK;
            v4df dr3_inv = dr2_inv * dr_inv;

            v4df m_j    = v4df(ep_j[j].mass);
            
#ifdef USE_QUAD
            v4df dr5_inv = dr3_inv * dr2_inv * v1p5;
            v4df dr7_inv = dr5_inv * dr2_inv;
            
            v4df qj_xx = v4df(ep_j[j].quad.xx);
            v4df qj_xy = v4df(ep_j[j].quad.xy);
            v4df qj_xz = v4df(ep_j[j].quad.xz);
            v4df qj_yy = v4df(ep_j[j].quad.yy);
            v4df qj_yz = v4df(ep_j[j].quad.yz);
            v4df qj_zz = v4df(ep_j[j].quad.zz);
            v4df qj_tr = v4df(ep_j[j].quad.getTrace());

            v4df qr_x = qj_xx*dpx_ij + qj_xy*dpy_ij + qj_xz*dpz_ij;
            v4df qr_y = qj_xy*dpx_ij + qj_yy*dpy_ij + qj_yz*dpz_ij;
            v4df qr_z = qj_xz*dpx_ij + qj_yz*dpy_ij + qj_zz*dpz_ij;

            v4df rqr = dpx_ij*qr_x + dpy_ij*qr_y + dpz_ij*qr_z;
            v4df A = m_j*dr3_inv - qj_tr*dr5_inv + v5*rqr*dr7_inv;
            v4df B = mv2*dr5_inv;

            pt_i -= m_j*dr_inv - v0p5*qj_tr*dr3_inv + rqr*dr5_inv;
            ax_i -= A*dpx_ij + B*qr_x;
            ay_i -= A*dpy_ij + B*qr_y;
            az_i -= A*dpz_ij + B*qr_z;
#else       
            v4df dg2_ij = m_j * dr3_inv;

            pt_i -= m_j    * dr_inv;
            ax_i -= dpx_ij * dg2_ij;
            ay_i -= dpy_ij * dg2_ij;
            az_i -= dpz_ij * dg2_ij;
#endif
        }

        PS::F64 buf_ax[nvector] __attribute__((aligned(32)));
        PS::F64 buf_ay[nvector] __attribute__((aligned(32)));
        PS::F64 buf_az[nvector] __attribute__((aligned(32)));
        PS::F64 buf_pt[nvector] __attribute__((aligned(32)));

        ax_i.store(buf_ax);
        ay_i.store(buf_ay);
        az_i.store(buf_az);
        pt_i.store(buf_pt);

        for(PS::S32 ii = 0; ii < nii; ii++) {
            force[i+ii].acc[0] += buf_ax[ii];
            force[i+ii].acc[1] += buf_ay[ii];
            force[i+ii].acc[2] += buf_az[ii];
            force[i+ii].phi    += buf_pt[ii];
        }
    }
}

template <class Tpsys>
void correctForceLong(Tpsys & pp,
#ifdef USE_QUAD
                      PS::TreeForForceLong<ForceGrav, EPGrav, EPGrav>::QuadrupoleWithScatterSearch & tree_grav,
#else
                      PS::TreeForForceLong<ForceGrav, EPGrav, EPGrav>::MonopoleWithScatterSearch & tree_grav,                
#endif                      
                      std::vector<std::vector<PS::S32> > & n_list,
                      PS::S32 & nei_dist,
                      PS::S32 & nei_tot_loc,
                      PS::S32 * ex_rank)
{
    PS::S32 n_loc = pp.getNumberOfParticleLocal();
    PS::F64 eps2 = EPGrav::eps2;
    PS::S32 j_id = 0;
    PS::S32 j_rank = 0;
    PS::S32 nei_loc = 0;
    clearNeighborList(n_list);
    n_list.resize(n_loc);
    
    PS::F64 r_IN = EPGrav::gamma*EPGrav::r_CUT;
    PS::F64 r_CUT_r_IN_inv = 1./((1.-EPGrav::gamma)*EPGrav::r_CUT);

    nei_dist = 0;
    nei_tot_loc = 0;

    for ( PS::S32 i=0; i<PS::Comm::getNumberOfProc(); i++ ) ex_rank[i] = 0;
    ex_rank[PS::Comm::getRank()] = 1;

    PS::S32 nei_tot_loc_tmp = 0;
    PS::S32 nei_dist_tmp = 0;
#pragma omp parallel for private(j_id, j_rank, nei_loc) reduction (+:nei_tot_loc_tmp, nei_dist_tmp)
    for(PS::S32 i = 0; i < n_loc; i++){
        pp[i].inDomain = true;
        EPGrav* next = NULL;//Pointer Of Neighbor List
        nei_loc = tree_grav.getNeighborListOneParticle(pp[i], next);
        pp[i].neighbor = 0;

        PS::F64vec posi = pp[i].getPos();
        PS::F64vec acci= 0.0;
        PS::F64 phii = 0.0;

        for(PS::S32 j = 0; j < nei_loc; j++){
            j_id = (next+j)->id;
            if ( pp[i].id == j_id ) continue;
            
            PS::F64vec posj   = (next+j)->getPos();
            PS::F64vec dr  = posj - posi;
            PS::F64 dr2 = dr * dr;
            if ( dr2 == 0.0 ) {
                if ( pp[i].id < 0 && j_id > -1 ) {
                    pp[i].id = j_id;
                    continue;
                } else if ( j_id < 0 ) continue;
            }
            assert( dr2 != 0.0 );
            dr2 += eps2;
            
#ifdef USE_INDIVIDUAL_RADII
            PS::F64 r_search = mutualRadius(pp[i].r_search, (next+j)->r_search);
            if ( dr2 < r_search * r_search ){
#else
                assert( dr2 < EPGrav::r_SEARCH * EPGrav::r_SEARCH );
#endif
                n_list[i].push_back(j_id);
                pp[i].neighbor ++;

                j_rank = (next+j)->myrank;
                if(j_rank != pp[i].myrank){
                    pp[i].inDomain = false;
                    ex_rank[j_rank] = 1;
                }
#ifdef USE_INDIVIDUAL_RADII
            }
#endif

            PS::F64 massj  = (next+j)->getCharge();

            PS::F64 rij   = sqrt(dr2);
            PS::F64 r2inv = 1.0 / dr2;
            PS::F64 rinv  = 1./rij;
            PS::F64 r3inv = rinv * r2inv;

#ifdef USE_INDIVIDUAL_RADII
            PS::F64 r_out = mutualRadius(pp[i].r_out, (next+j)->r_out);
#else
            PS::F64 r_out = EPGrav::r_CUT;
#endif      
            PS::F64 W  = cutoff_W(rij, r_out);
            PS::F64 K  = cutoff_K(rij, r_out);
            PS::F64 K0 = std::min(std::max((rij - r_IN)*r_CUT_r_IN_inv, 0.), 1. );

            phii   -= massj * rinv * (W - K0);
            acci  += massj * r3inv * (K - K0) * dr;
        }
        
        nei_tot_loc_tmp += pp[i].neighbor;
        nei_dist_tmp ++;
            
        pp[i].acc  += acci;
        pp[i].phi  += phii;
    }
    
    nei_tot_loc += nei_tot_loc_tmp;
    nei_dist += nei_dist_tmp;
}

template <class Tpsys>
void correctForceLongInitial(Tpsys & pp,
#ifdef USE_QUAD
                             PS::TreeForForceLong<ForceGrav, EPGrav, EPGrav>::QuadrupoleWithScatterSearch & tree_grav,
#else
                             PS::TreeForForceLong<ForceGrav, EPGrav, EPGrav>::MonopoleWithScatterSearch & tree_grav,                
#endif
                             std::vector<std::vector<PS::S32> > & n_list,
                             PS::S32 & nei_dist,
                             PS::S32 & nei_tot_loc,
                             PS::S32 * ex_rank)
{
    PS::S32 n_loc = pp.getNumberOfParticleLocal();
    PS::F64 eps2 = EPGrav::eps2;
    PS::S32 j_id = 0;
    PS::S32 j_rank = 0;
    PS::S32 nei_loc = 0;
    clearNeighborList(n_list);
    n_list.resize(n_loc);
    
    PS::F64 r_IN = EPGrav::gamma*EPGrav::r_CUT;
    PS::F64 r_CUT_r_IN_inv = 1./((1.-EPGrav::gamma)*EPGrav::r_CUT);

    nei_dist = 0;
    nei_tot_loc = 0;

    for ( PS::S32 i=0; i<PS::Comm::getNumberOfProc(); i++ ) ex_rank[i] = 0;
    ex_rank[PS::Comm::getRank()] = 1;

    PS::S32 nei_tot_loc_tmp = 0;
    PS::S32 nei_dist_tmp = 0;
#pragma omp parallel for private(j_id, j_rank, nei_loc) reduction (+:nei_tot_loc_tmp, nei_dist_tmp)
    for(PS::S32 i = 0; i < n_loc; i++){
        pp[i].acc_d = 0.0;
        pp[i].jerk = 0.0;
        pp[i].phi_d = 0.0;
        pp[i].phi_s = 0.0;

        calcStarGravity(pp[i]);
        
        pp[i].inDomain = true;
        EPGrav* next = NULL;//Pointer Of Neighbor List
        nei_loc = tree_grav.getNeighborListOneParticle(pp[i], next);
        pp[i].neighbor = 0;

        PS::F64vec posi = pp[i].getPos();
        PS::F64vec veli = pp[i].vel;
        PS::F64vec acci = 0.0;
        PS::F64vec acc_di = 0.0;
        PS::F64vec jerki = 0.0;
        PS::F64 phii = 0.0;
        PS::F64 phi_di = 0.0;

        for(PS::S32 j = 0; j < nei_loc; j++){
            j_id = (next+j)->id;
            if ( pp[i].id == j_id ) continue;
            
            PS::F64vec posj   = (next+j)->getPos();
            PS::F64vec dr  = posj - posi;
            PS::F64 dr2 = dr * dr;
            if ( dr2 == 0.0 ) {
                if ( pp[i].id < 0 && j_id > -1 ) {
                    pp[i].id = j_id;
                    continue;
                } else if ( j_id < 0 ) continue;
            }
            assert( dr2 != 0.0 );
            dr2 += eps2;

#ifdef USE_INDIVIDUAL_RADII
            PS::F64 r_search = mutualRadius(pp[i].r_search, (next+j)->r_search);
            if ( dr2 < r_search * r_search ){
#else
                assert( dr2 < EPGrav::r_SEARCH * EPGrav::r_SEARCH );
#endif
                n_list[i].push_back(j_id);
                pp[i].neighbor ++;

                j_rank = (next+j)->myrank;
                if(j_rank != pp[i].myrank){
                    pp[i].inDomain = false;
                    ex_rank[j_rank] = 1;
                }
#ifdef USE_INDIVIDUAL_RADII
            }
#endif
            PS::F64vec velj   = (next+j)->vel;
            PS::F64vec dv  = velj - veli;
            PS::F64 drdv  = dr*dv;
            PS::F64 massj  = (next+j)->getCharge();

            PS::F64 rij   = sqrt(dr2);
            PS::F64 r2inv = 1./dr2;
            PS::F64 rinv  = 1./rij;
            PS::F64 r3inv = rinv * r2inv;
            PS::F64 r5inv = r3inv * r2inv;

#ifdef USE_INDIVIDUAL_RADII
            PS::F64 r_out = mutualRadius(pp[i].r_out, (next+j)->r_out);
#else
            PS::F64 r_out = EPGrav::r_CUT;
#endif
            PS::F64 W  = cutoff_W(rij, r_out);
            PS::F64 K  = cutoff_K(rij, r_out);
            PS::F64 dK = cutoff_dK(rij, drdv, r_out);
            PS::F64 K0 = std::min(std::max((rij - r_IN)*r_CUT_r_IN_inv, 0.), 1. );
            
            phii   -= massj * rinv * (W - K0);
            phi_di -= massj * rinv * (1.-W);
            acci   += massj * r3inv * (K - K0) * dr;
            acc_di += massj * r3inv * dr * (1.-K);
            jerki  += massj * ( (r3inv*dv -3.*r5inv*(drdv)*dr)*(1.-K) - r3inv*dr*dK );
        }
        
        nei_tot_loc_tmp += pp[i].neighbor;
        nei_dist_tmp ++;
        
        pp[i].acc  += acci;
        pp[i].phi  += phii;
        pp[i].acc_d += acc_di;
        pp[i].phi_d += phi_di;
        pp[i].jerk += jerki;

        pp[i].calcDeltatInitial();
    }

    nei_tot_loc += nei_tot_loc_tmp;
    nei_dist += nei_dist_tmp;
}


