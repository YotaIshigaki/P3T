#pragma once

class HardSystem{
public:
    std::vector<PS::S32> list_iso;
    std::vector<std::vector<std::pair<bool,PS::S32> > > list_multi;
    std::vector<std::vector<FPHard> > ptcl_multi;
    std::map<PS::S32,PS::S32> mp_cluster;  // cluster ID -> cluster adress in HardSystem

    std::vector<Collision> collision_list;
    std::vector<std::pair<PS::S32,PS::S32> > frag_list;

    PS::S32 n_col;
    PS::S32 n_frag;
    PS::F64 edisp;
    PS::F64 edisp_d;

    static PS::F64 f;

    PS::S32 getNumberOfClusterLocal() const { return ptcl_multi.size(); }
    PS::S32 getNumberOfClusterGlobal() const {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        PS::S32 n_cluster_loc = ptcl_multi.size();
        return PS::Comm::getSum(n_cluster_loc);
#else
        return ptcl_multi.size();
#endif
    }
    PS::S32 getNumberOfIsolatedParticleLocal() const { return list_iso.size(); }
    PS::S32 getNumberOfIsolatedParticleGlobal() const {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        PS::S32 n_ptcl_loc = list_iso.size();
        return PS::Comm::getSum(n_ptcl_loc);
#else
        return list_iso.size();
#endif
    }
    PS::S32 getNumberOfParticleLocal(){
        PS::S32 n = 0;
        for ( PS::S32 i=0; i<ptcl_multi.size(); i++ ) n += ptcl_multi[i].size();
        n += list_iso.size();
        return n;
    }
    PS::S32 getNumberOfParticleGlobal() {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        PS::S32 n_ptcl_loc = getNumberOfParticleLocal();
        return PS::Comm::getSum(n_ptcl_loc);
#else
        return getNumberOfParticleLocal();
#endif
    }
    PS::S32 getNumberOfParticleInLargestClusterLocal(){
        PS::S32 n = 0;
        for ( PS::S32 i=0; i<ptcl_multi.size(); i++ )
            if ( n < ptcl_multi[i].size() ) n = ptcl_multi[i].size();
        return n;
    }
    PS::S32 getNumberOfParticleInLargestClusterGlobal() {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        PS::S32 n_loc = getNumberOfParticleInLargestClusterLocal();
        return PS::Comm::getMaxValue(n_loc);
#else
        return getNumberOfParticleInLargestClusterLocal();
#endif
    }
    
    PS::S32 getNumberOfFragmentLocal() const { return n_frag; }
    PS::S32 getNumberOfFragmentGlobal() const {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        PS::S32 n_frag_ = n_frag;
        return PS::Comm::getSum(n_frag_);
#else
        return n_frag;
#endif
    }
    PS::S32 getNumberOfCollisionLocal() const { return n_col; }
    PS::S32 getNumberOfCollisionGlobal() const {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        PS::S32 n_col_ = n_col;
        return PS::Comm::getSum(n_col_);
#else
        return n_col;
#endif
    }
    PS::F64 getEnergyDissipationLocal() const { return edisp; }
    PS::F64 getEnergyDissipationGlobal() const {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        PS::F64 edisp_ = edisp;
        return PS::Comm::getSum(edisp_);
#else
        return edisp;
#endif
    }
    PS::F64 getHardEnergyDissipationLocal() const { return edisp_d; }
    PS::F64 getHardEnergyDissipationGlobal() const {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        PS::F64 edisp_d_ = edisp_d;
        return PS::Comm::getSum(edisp_d_);
#else
        return edisp_d;
#endif
    }
    
    void showParticleID() const {
        for ( PS::S32 i=0; i<ptcl_multi.size(); i++ ){
            bool flag = ptcl_multi.at(i).at(0).inDomain;
            std::cout << "Rank " << PS::Comm::getRank()
                      << " Cluster " << ptcl_multi.at(i).at(0).id_cluster
                      << " (" << flag << "): ";
            for ( PS::S32 j=0; j<ptcl_multi.at(i).size(); j++ ){
                std::cout << " " << ptcl_multi.at(i).at(j).id;
                assert( flag == ptcl_multi.at(i).at(j).inDomain );
            }
            std::cout << std::endl;
        }
    }

    void clear(){
        list_iso.clear();
        list_multi.clear();        
        ptcl_multi.clear();
        mp_cluster.clear();
        collision_list.clear();
        frag_list.clear();

        n_col = n_frag = 0;
        edisp = edisp_d = 0.;
    }

    template <class Tpsys, class Tp>
    PS::S32 makeList(Tpsys & pp,
                     PS::S32 n_pp,
                     Tp * & ex_pp,
                     PS::S32 n_ex_pp);

    template <class Tpsys, class Tp>
    PS::S32 timeIntegrate(Tpsys & pp,
                          Tp * & ex_pp,
                          std::vector<std::vector<PS::S32> > & n_list,
                          PS::S32 * ex_nei_list,
                          PS::S32 * ex_nei_adr_list,
                          const PS::S32 istep);

    static void rewriteFragmentID(PS::S32 * & id_frag_list,
                                  Collision * & col_list,
                                  PS::S32 n_col_tot,
                                  PS::S32 n_frag_tot,
                                  PS::S32 & id_next);
    template <class Tpsys>
    PS::S32 addFragment2ParticleSystem(Tpsys & pp,
                                       PS::S32 & id_next,
                                       std::ofstream & fp); 
};

PS::F64 HardSystem::f = 1.;


template <class Tpsys, class Tp>
inline PS::S32 HardSystem::makeList(Tpsys & pp,
                                    PS::S32 n_pp,
                                    Tp * & ex_pp,
                                    PS::S32 n_ex_pp)
{   
    PS::S32 n_ptcl_loc = 0;
    PS::S32 n_cluster_loc = 0;
    PS::S32 tmp = 0;
        
    for ( PS::S32 i=0; i<n_pp; i++ ){
        if ( pp[i].neighbor ) {
            if ( !pp[i].inDomain && n_ex_pp == 0 ) continue;
            auto itr = mp_cluster.find(pp[i].id_cluster);                    
            if ( itr == mp_cluster.end() ){
                mp_cluster[pp[i].id_cluster] = tmp;
                std::vector<std::pair<bool,PS::S32> > tmp_v{std::make_pair(true, i)};
                list_multi.push_back(tmp_v);
                tmp ++;
                assert( list_multi.size() == tmp );
            } else {
                assert( mp_cluster.at(pp[i].id_cluster) == itr->second );
                list_multi.at(itr->second).push_back(std::make_pair(true, i));
            }
            n_ptcl_loc ++;
        } else {
            list_iso.push_back(i);
            n_ptcl_loc ++;
        }
    }
    
    //ex_nei_adr_list[0] = 0;
    for ( PS::S32 i=0; i<n_ex_pp; i++ ){
        assert( ex_pp[i].neighbor > 0 );
        assert( ex_pp[i].inDomain == false );
        auto itr = mp_cluster.find(ex_pp[i].id_cluster);                    
        if ( itr == mp_cluster.end() ){
            mp_cluster[ex_pp[i].id_cluster] = tmp;
            std::vector<std::pair<bool,PS::S32> > tmp_v{std::make_pair(false, i)};
            list_multi.push_back(tmp_v);
            tmp ++;
            assert( list_multi.size() == tmp );
        } else {
            assert( mp_cluster.at(ex_pp[i].id_cluster) == itr->second );
            list_multi.at(itr->second).push_back(std::make_pair(false, i));
        }
        n_ptcl_loc ++;
        //if ( i>0 )  ex_nei_adr_list[i] = ex_nei_adr_list[i-1] + ex_pp[i-1].neighbor;
    }
    n_cluster_loc = tmp;
    assert( list_multi.size() == tmp );
    //assert( ex_nei_adr_list[n_ex_pp-1]+ex_pp[n_ex_pp-1].neighbor == n_ex_nei );
    ptcl_multi.resize(n_cluster_loc);

    return n_ptcl_loc;
}

template <class Tpsys, class Tp>
inline PS::S32 HardSystem::timeIntegrate(Tpsys & pp,
                                         Tp * & ex_pp,
                                         std::vector<std::vector<PS::S32> > & n_list,
                                         PS::S32 * ex_nei_list,
                                         PS::S32 * ex_nei_adr_list,
                                         const PS::S32 istep)
{
    PS::S32 n_all = list_multi.size() + list_iso.size();
    PS::S32 n_ptcl_loc = 0;
    
#pragma omp parallel for reduction(+:n_ptcl_loc) schedule (dynamic)
    for ( PS::S32 ii=0; ii<n_all; ii++ ){
        if ( ii<list_multi.size() ){
            PS::S32 i = ii;
            //for ( PS::S32 i=0; i<list_multi.size(); i++ ){
            std::map<PS::S32,PS::S32> id2adr;
            PS::S32 n_p = list_multi.at(i).size();
            PS::S32 id_cluster = 0;
            ptcl_multi[i].clear();
            id2adr.clear();

            // Add Particle To Hard System
            ptcl_multi[i].reserve(n_p);
            for ( PS::S32 j=0; j<n_p; j++ ){
                std::pair<bool,PS::S32> adr = list_multi.at(i).at(j);
                if ( adr.first ) {
                    ptcl_multi[i].push_back(FPHard(pp[adr.second]));
                    ptcl_multi[i][j].copyList(n_list.at(adr.second));
                } else {
                    ptcl_multi[i].push_back(FPHard(ex_pp[adr.second]));
                    ptcl_multi[i][j].copyList(&ex_nei_list[ex_nei_adr_list[adr.second]]);
                }
                if ( j==0 ) id_cluster = ptcl_multi[i][j].id_cluster;
                assert ( ptcl_multi[i][j].id_cluster == id_cluster );
                id2adr[ptcl_multi[i][j].id] = j;
            }

            // Make Neighbor List
            for ( PS::S32 j=0; j<n_p; j++ ){
                ptcl_multi[i][j].n_hard_list.clear();
                ptcl_multi[i][j].n_hard_list.reserve(ptcl_multi[i][j].neighbor);
                for ( PS::S32 k=0; k<ptcl_multi[i][j].neighbor; k++){
                    ptcl_multi[i][j].n_hard_list.push_back(id2adr.at(ptcl_multi[i][j].n_list.at(k)));
                }
            }
        
            PS::S32 n_col_tmp  = 0;
            PS::S32 n_frag_tmp = 0;
            PS::F64 edisp_tmp   = 0.;
            PS::F64 edisp_d_tmp = 0.;
            std::vector<Collision> collision_list_tmp;
            timeIntegrate_multi(ptcl_multi[i], (istep-1)*FPGrav::dt_tree, istep*FPGrav::dt_tree, f,
                                n_col_tmp, n_frag_tmp, edisp_tmp, edisp_d_tmp, collision_list_tmp);
            if ( n_col_tmp > 0 ){
#pragma omp critical
                {
                    n_col  += n_col_tmp;
                    n_frag += n_frag_tmp;
                    edisp   += edisp_tmp;
                    edisp_d += edisp_d_tmp;
                    for ( PS::S32 j=0; j<n_col_tmp; j++ )
                        collision_list.push_back(collision_list_tmp.at(j));
                    for ( PS::S32 j=0; j<n_frag_tmp; j++ )
                        frag_list.push_back(std::make_pair(i, ptcl_multi[i].size()-n_frag_tmp+j));
                }
            }
                
            for ( PS::S32 j=0; j<ptcl_multi[i].size()-n_frag_tmp; j++ ){
                std::pair<bool,PS::S32> adr = list_multi.at(i).at(j);
                if ( adr.first ) {
                    if ( !ptcl_multi[i][j].isDead ) assert ( pp[adr.second].id == ptcl_multi[i][j].id );
                    pp[adr.second] = FPGrav(ptcl_multi[i][j]);
                } else {
                    if ( !ptcl_multi[i][j].isDead ) assert ( ex_pp[adr.second].id == ptcl_multi[i][j].id );
                    ex_pp[adr.second] = FPGrav(ptcl_multi[i][j]);
                }
                n_ptcl_loc ++;
            }
            
        } else {
            PS::S32 i = ii - list_multi.size();
            //#pragma omp parallel for reduction(+:n_ptcl_loc)
            //for ( PS::S32 i=0; i<list_iso.size(); i++ ){
            if ( pp[list_iso[i]].getSemimajorAxis() > 0 && FPGrav::eps2 == 0. ){
                //if (0) {
                timeIntegrateKepler_isolated(pp[list_iso[i]], (istep-1)*FPGrav::dt_tree, istep*FPGrav::dt_tree);
            } else {
                FPHard pi = FPHard(pp[list_iso[i]]);
                timeIntegrate_isolated(pi, (istep-1)*FPGrav::dt_tree, istep*FPGrav::dt_tree);
                pp[list_iso[i]] = FPGrav(pi);
            }
            n_ptcl_loc ++;
        }
    }
    
    return n_ptcl_loc;
}

inline void HardSystem::rewriteFragmentID(PS::S32 * & id_frag_list,
                                          Collision * & col_list,
                                          PS::S32 n_col_tot,
                                          PS::S32 n_frag_tot,
                                          PS::S32 & id_next)
{
    std::map<PS::S32, PS::S32> id_old2new;
    id_old2new.clear();
        
    for ( PS::S32 i=0; i<n_col_tot; i++ ){
        PS::S32 n_fragi  = col_list[i].getNumberOfFragment();
        PS::S32 id_fragi = col_list[i].getFragmentID();
        if ( n_fragi == 0 ) continue;
        for ( PS::S32 j=0; j<n_fragi; j++ ){
            id_old2new[id_fragi - j] = id_next;
            id_next ++;
        }
    }
    
    for ( PS::S32 i=0; i<n_frag_tot; i++ )
        if ( id_frag_list[i] < 0 ) id_frag_list[i] = id_old2new.at(id_frag_list[i]);
    for ( PS::S32 i=0; i<n_col_tot; i++ )
        col_list[i].setNewFragmentID(id_old2new);
}
template <class Tpsys>
inline PS::S32 HardSystem::addFragment2ParticleSystem(Tpsys & pp,
                                                      PS::S32 & id_next,
                                                      std::ofstream & fp)
{
    PS::S32 * n_col_list     = nullptr;
    PS::S32 * n_frag_list    = nullptr;
    Collision * col_list_tot = nullptr;
    PS::S32 * id_frag_list   = nullptr;
    PS::S32 * id_frag_loc    = nullptr;
    PS::S32 * col_recv  = nullptr;
    PS::S32 * frag_recv = nullptr;

    id_frag_loc = new PS::S32[n_frag];
    for ( PS::S32 i=0; i<n_frag; i++ ){
        std::pair<PS::S32, PS::S32> adr = frag_list.at(i);
        id_frag_loc[i] = ptcl_multi[adr.first][adr.second].id;
    }
    if ( PS::Comm::getRank() == 0 ){
        n_col_list  = new PS::S32[PS::Comm::getNumberOfProc()];
        n_frag_list = new PS::S32[PS::Comm::getNumberOfProc()];
        col_recv    = new PS::S32[PS::Comm::getNumberOfProc()];
        frag_recv   = new PS::S32[PS::Comm::getNumberOfProc()];
        col_recv[0]  = 0;
        frag_recv[0] = 0;
    }
    // Send Number of Collision & Fragments
    PS::Comm::gather(&n_col, 1, n_col_list);
    PS::Comm::gather(&n_frag, 1, n_frag_list);

    PS::S32 n_col_tot = 0;
    PS::S32 n_frag_tot = 0;
    if ( PS::Comm::getRank() == 0 ){
        PS::S32 tmp_col = 0;
        PS::S32 tmp_frag = 0;
        for ( PS::S32 i=1; i<PS::Comm::getNumberOfProc(); i++ ){
            tmp_col  += n_col_list[i-1];
            tmp_frag += n_frag_list[i-1];
            col_recv[i]  = tmp_col;
            frag_recv[i] = tmp_frag;
        }
        tmp_col  += n_col_list[PS::Comm::getNumberOfProc()-1];
        tmp_frag += n_frag_list[PS::Comm::getNumberOfProc()-1];
        col_list_tot = new Collision[tmp_col];
        id_frag_list = new PS::S32[tmp_frag];
        n_col_tot = tmp_col;
        n_frag_tot = tmp_frag;
    }
    // Send Collision Information & Fragments ID
    PS::Comm::gatherV(&collision_list[0], n_col,  col_list_tot, n_col_list,  col_recv);
    PS::Comm::gatherV(id_frag_loc,        n_frag, id_frag_list, n_frag_list, frag_recv);

    // Rewrite Fragments ID
    if ( PS::Comm::getRank() == 0 ){
        rewriteFragmentID(id_frag_list, col_list_tot, n_col_tot, n_frag_tot, id_next);
    }

    // Return Fragments ID
    PS::Comm::scatterV(id_frag_list, n_frag_list, frag_recv, id_frag_loc, n_frag);
    assert( frag_list.size() == n_frag );
    for ( PS::S32 i=0; i<n_frag; i++ ){
        std::pair<PS::S32, PS::S32> adr = frag_list.at(i);
        ptcl_multi[adr.first][adr.second].id = id_frag_loc[i];
        pp.addOneParticle(ptcl_multi[adr.first][adr.second]);
    }
    PS::Comm::broadcast(&id_next, 1);

    if ( PS::Comm::getRank() == 0 ){
        // Output Collision Information
        for ( PS::S32 i=0; i<n_col_tot; i++ ) col_list_tot[i].write2File(fp);
            
        delete [] n_col_list;
        delete [] n_frag_list;
        delete [] col_list_tot;
        delete [] id_frag_list;
        delete [] col_recv;
        delete [] frag_recv;
    }
    delete [] id_frag_loc;

    return n_frag_tot;
}





template <class Tpsys>
void createNeighborCluster(Tpsys & pp,
                           std::map<PS::S32, PS::S32> & id2id_loc,
                           std::vector<std::vector<PS::S32> > n_list)
{
    const PS::S32 n_loc = pp.getNumberOfParticleLocal();
    PS::S32 nei = 0;
    PS::S32 j_id = 0;
    PS::S32 id_cluster = 0;
    //List Of Particles Which Have Neighbors Out Of Domain
    std::vector<PS::S32> out_of_domain_list;
    assert( n_list.size() == n_loc );
    
    for(PS::S32 i=0; i<n_loc; i++){
        assert( pp[i].neighbor == n_list[i].size() );
        nei = pp[i].neighbor;
        id_cluster = pp[i].id;
        pp[i].id_cluster = pp[i].id;
        id2id_loc[pp[i].id] = i;

        if(nei == 0) continue;
        for(PS::S32 j=0; j<nei; j++){                
            j_id = n_list[i].at(j);
            if( id_cluster > j_id ) id_cluster = j_id;
        }
        pp[i].id_cluster = id_cluster;
        if( pp[i].inDomain == false ) out_of_domain_list.push_back(pp[i].id_local);
    }

    PS::S32 j_id_cluster = 0;
    bool check = true;
    while( check ){
        check = false;
        for(PS::S32 i=0; i<n_loc; i++){
            nei = pp[i].neighbor;
            id_cluster = pp[i].id_cluster;

            if(nei == 0) continue;
            for(PS::S32 j=0; j<nei; j++){
                auto itr = id2id_loc.find(n_list[i].at(j));
                if ( itr == id2id_loc.end() ) continue;
                j_id = itr->second;
                j_id_cluster = pp[j_id].id_cluster;
                if( id_cluster > j_id_cluster ) id_cluster = j_id_cluster;
            }
            if( pp[i].id_cluster != id_cluster ) check = true;
            pp[i].id_cluster = id_cluster;
            assert( pp[i].id >= id_cluster );
        }
    }

    if( out_of_domain_list.size() != 0 ){
        PS::S32 out_cluster_id = 0;
        PS::S32 out_particle_local = 0;
        nei = out_of_domain_list.size();
        for(PS::S32 i=0; i<n_loc; i++){
            id_cluster = pp[i].id_cluster;
            for(PS::S32 j=0; j<nei; j++){
                out_particle_local = out_of_domain_list.at(j);
                out_cluster_id = pp[out_particle_local].id_cluster;
                if( id_cluster == out_cluster_id ){
                    pp[i].inDomain = false;
                }
            }
        }
    }
}
template <class Tp, class Tpsys>
void createNeighborCluster_OutOfDomain(Tpsys & pp,
                                       PS::S32 n_pp,
                                       Tp * & ex_pp,
                                       PS::S32 n_ex_pp,
                                       std::map<PS::S32, PS::S32> & id2id_loc,
                                       std::vector<std::vector<PS::S32> > & n_list,
                                       PS::S32 * & ex_nei_list,
                                       PS::S32 * & ex_nei_adr_list,
                                       //PS::S32 n_ex_nei,
                                       PS::S32 recv_rank)
{
    std::vector<PS::S32> out_of_domain_list;
    std::map<PS::S32, PS::S32>  id2id_ex;
    out_of_domain_list.clear();
    id2id_ex.clear();
    for (PS::S32 i=0; i<n_pp; i++){
        if ( !pp[i].inDomain ) out_of_domain_list.push_back(i);
    }
    for (PS::S32 i=0; i<n_ex_pp; i++) id2id_ex[ex_pp[i].id] = i;
    
    PS::S32 id_cluster = 0;
    PS::S32 nei = 0;
    PS::S32 j_id = 0;
    PS::S32 j_id_cluster = 0;
    bool check = true;
    while( check ){
        check = false;
        for(PS::S32 i=0; i<out_of_domain_list.size(); i++){
            PS::S32 i_id = out_of_domain_list.at(i);
            nei = pp[i_id].neighbor;
            id_cluster = pp[i_id].id_cluster;
            assert (nei > 0);
            
            for(PS::S32 j=0; j<nei; j++){
                PS::S32 id_nei = n_list[i_id].at(j);
                auto itr = id2id_loc.find(id_nei);
                auto itr_ex = id2id_ex.find(id_nei);
                if ( itr != id2id_loc.end() ) {
                    assert ( itr_ex == id2id_ex.end() );
                    j_id = itr->second;
                    j_id_cluster = pp[j_id].id_cluster;
                    if( id_cluster > j_id_cluster ) id_cluster = j_id_cluster;
                    assert ( pp[j_id].inDomain == false );
                } else {
                    assert ( itr_ex != id2id_ex.end() );
                    j_id = itr_ex->second;
                    j_id_cluster = ex_pp[j_id].id_cluster;
                    if( id_cluster > j_id_cluster ) id_cluster = j_id_cluster;
                    assert ( ex_pp[j_id].inDomain == false );
                }
            }
            if( pp[i_id].id_cluster != id_cluster ) check = true;
            pp[i_id].id_cluster = id_cluster;
        }
        //PS::S32 nei_tmp = 0;
        for(PS::S32 i=0; i<n_ex_pp; i++){
            nei = ex_pp[i].neighbor;
            id_cluster = ex_pp[i].id_cluster;
            assert (nei > 0);
            
            for(PS::S32 j=0; j<nei; j++){
                PS::S32 id_nei = ex_nei_list[ex_nei_adr_list[i]+j];
                auto itr = id2id_loc.find(id_nei);
                auto itr_ex = id2id_ex.find(id_nei);
                if ( itr != id2id_loc.end() ) {
                    j_id = itr->second;
                    j_id_cluster = pp[j_id].id_cluster;
                    if( id_cluster > j_id_cluster ) id_cluster = j_id_cluster;
                    assert ( pp[j_id].inDomain == false );
                } else {
                    assert ( itr_ex != id2id_ex.end() );
                    j_id = itr_ex->second;
                    j_id_cluster = ex_pp[j_id].id_cluster;
                    if( id_cluster > j_id_cluster ) id_cluster = j_id_cluster;
                    assert ( ex_pp[j_id].inDomain == false );
                }
            }
            if( ex_pp[i].id_cluster != id_cluster ) check = true;
            ex_pp[i].id_cluster = id_cluster;
            assert( ex_pp[i].id >= id_cluster );

            //nei_tmp += nei;
        }
        //assert( nei_tmp == n_ex_nei );
    }
}


