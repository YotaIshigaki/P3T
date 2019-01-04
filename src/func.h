#pragma once

template <class Tpsys>
void calcMeanMass(Tpsys & pp,
                  PS::F64 & m_mean,
                  PS::F64 & m_max)
{
    const PS::S32 n_loc = pp.getNumberOfParticleLocal();
    const PS::S32 n_glb = pp.getNumberOfParticleGlobal();
    PS::F64 m_sum_loc = 0.;
    PS::F64 m_max_loc = 0.;
    
    for (PS::S32 i=0; i<n_loc; i++ ){
        m_sum_loc += pp[i].mass;
        if ( pp[i].mass > m_max_loc ) m_max_loc = pp[i].mass;
    }
    m_mean = PS::Comm::getSum(m_sum_loc) / n_glb;
    m_max = PS::Comm::getMaxValue(m_max_loc);
}

template <class Tpsys>
void makeSnap(Tpsys & pp,
              PS::F64 time_sys,
              Energy e_init,
              Energy e_now,
              const char * dir_name,
              const PS::S32 isnap)
{
    FileHeader header(pp.getNumberOfParticleGlobal(), time_sys, e_init, e_now);
    char filename[256];
    sprintf(filename, "%s/snap%06d.dat", dir_name, isnap);
    pp.writeParticleAscii(filename, header);
}

template <class Tpsys>
void outputStep(Tpsys & pp,
                PS::F64 time_sys,
                Energy e_init,
                Energy e_now,
                PS::F64 de,
                PS::S32 n_col_tot,
                PS::S32 n_frag_tot,
                const char * dir_name,
                const PS::S32 isnap,
                std::ofstream & fout_eng,
                Wtime wtime,
                PS::S32 n_largestcluster,
                PS::S32 n_cluster,
                PS::S32 n_isoparticle,
                bool bSnap=true)
{
    const PS::S32 n_tot = pp.getNumberOfParticleGlobal();

    if ( bSnap ) makeSnap(pp, time_sys, e_init, e_now, dir_name, isnap);

#ifdef OUTPUT_DETAIL
    PS::F64 m_mean = 0.;
    PS::F64 m_max = 0.;
    calcMeanMass(pp, m_mean, m_max);
#endif

    if(PS::Comm::getRank() == 0 && bSnap){
        //PS::F64 de =  e_now.calcEnergyError(e_init);
        //PS::F64 de_tmp = sqrt(de*de);
        //if( de_tmp > de_max ) de_max = de_tmp;
        fout_eng  << std::fixed<<std::setprecision(8)
                  << time_sys << "\t" << n_tot << "\t"
                  << std::scientific<<std::setprecision(15)
                  << e_now.etot << "\t" << de << "\t"
                  << n_largestcluster << "\t" << n_cluster << "\t" << n_isoparticle
#ifdef OUTPUT_DETAIL
                  << "\t" << m_max << "\t" << m_mean 
#endif
#ifdef CALC_WTIME
                  << "\t" << wtime.soft_step << "\t" << wtime.hard_step << "\t"
                  << wtime.calc_soft_force_step << "\t" << wtime.neighbor_search_step << "\t"
                  << wtime.calc_hard_force_step << "\t"
                  << wtime.create_cluster_step << "\t" << wtime.communication_step << "\t"
                  << wtime.output_step 
#endif
                  <<std::endl;
    }
}

template <class Tpsys>
void inputIDLocalAndMyrank(Tpsys & pp)
{
    const PS::S32 n_loc = pp.getNumberOfParticleLocal();
    PS::S32 myrank = PS::Comm::getRank();
#pragma omp parallel for
    for(PS::S32 i=0; i<n_loc; i++){
        pp[i].id_local = i;
        pp[i].myrank = myrank;
    }
}

void decideExchangeRank(PS::S32 * ex_rank,
                        PS::S32 * ex_rank_tmp,
                        PS::S32 & n_ex_rank,
                        PS::S32 & recv_rank)
{
    const PS::S32 myrank = PS::Comm::getRank();
    const PS::S32 n_Proc = PS::Comm::getNumberOfProc();
    bool check;
    bool check_tot;
    PS::S32 loop = 0;

    std::vector<PS::S32> connected_rank_list;
    connected_rank_list.clear();
    for ( PS::S32 i=0; i<n_Proc; i++ ){
        if ( ex_rank[i] ) connected_rank_list.push_back(i);
    }

    //PS::Comm::barrier();
    do {
        PS::S32 tmp = 0;
        for ( PS::S32 i=0; i<n_Proc; i++ ){
            if ( ex_rank[i] ) {
                if ( tmp == 0 ) recv_rank = i;
                tmp ++;
            }
        }
        n_ex_rank = tmp;
        check = true;
       
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        for ( PS::S32 i=0; i<connected_rank_list.size(); i++ ){
            PS::S32 rank = connected_rank_list.at(i);
            if ( rank == myrank ) continue;
            PS::S32 TAG = loop + 10;
            MPI_Status status;
           
            //MPI::COMM_WORLD.Sendrecv(ex_rank,   n_Proc, MPI::BOOL, rank, TAG,
            //                         ex_rank_tmp, n_Proc, MPI::BOOL, rank, TAG);
            MPI_Sendrecv(ex_rank,     n_Proc, MPI_INT, rank, TAG,
                         ex_rank_tmp, n_Proc, MPI_INT, rank, TAG, MPI_COMM_WORLD, &status);
            
            for ( PS::S32 j=0; j<n_Proc; j++ ){
                //if ( ( ex_rank[j] && !ex_rank_tmp[j] ) || ( !ex_rank[j] && ex_rank_tmp[j] ) ){
                if ( ex_rank[j] != ex_rank_tmp[j] ){
                    check = false;
                    ex_rank[j] = 1;
                }
            }
        }
        check_tot = PS::Comm::synchronizeConditionalBranchAND(check);
#else
        assert ( connected_rank_list.size() == 1 );
        check_tot = true;
#endif
        loop ++;
    } while ( !check_tot );
}

template <class Tpsys>
void countSendParticle(Tpsys & pp,
                       PS::S32 & ex_tot,
                       PS::S32 & ex_loc,
                       PS::S32 & ex_nei_loc)
{
    const PS::S32 n_loc = pp.getNumberOfParticleLocal();
    ex_tot = 0;
    ex_loc = 0;
    ex_nei_loc = 0;
    PS::S32 ex_loc_tmp = 0;
    PS::S32 ex_nei_loc_tmp = 0;
#pragma omp parallel for reduction (+:ex_loc_tmp, ex_nei_loc_tmp)
    for(PS::S32 i=0; i<n_loc; i++){
        if(pp[i].inDomain == false){
            ex_loc_tmp ++;
            ex_nei_loc_tmp += pp[i].neighbor;
        }
    }
    ex_loc += ex_loc_tmp;
    ex_nei_loc += ex_nei_loc_tmp;
    ex_tot = PS::Comm::getSum(ex_loc);
}

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
void SendRecvOutOfDomainParticleNumber(PS::S32 * ex_rank,
                                       PS::S32 recv_rank,
                                       PS::S32 ex_loc,
                                       PS::S32 ex_nei_loc,
                                       PS::S32 * & n_ex_ptcl_,
                                       PS::S32 * & n_ex_nei_)
{
    PS::S32 TAG = 0;
    if ( ex_loc ) {
        if ( PS::Comm::getRank() != recv_rank ){
            assert( PS::Comm::getRank() > recv_rank );
            //MPI::COMM_WORLD.Send(&ex_loc,     1, PS::GetDataType(ex_loc),     recv_rank, TAG);
            //MPI::COMM_WORLD.Send(&ex_nei_loc, 1, PS::GetDataType(ex_nei_loc), recv_rank, TAG+1);
            MPI_Send(&ex_loc,     1, PS::GetDataType(ex_loc),     recv_rank, TAG,   MPI_COMM_WORLD);
            MPI_Send(&ex_nei_loc, 1, PS::GetDataType(ex_nei_loc), recv_rank, TAG+1, MPI_COMM_WORLD);
        } else {               
            PS::S32 tmp = 0;
            for(PS::S32 i=0; i<PS::Comm::getNumberOfProc(); ++i){
                if ( ex_rank[i] ){
                    if ( i == PS::Comm::getRank() ) continue;
                    MPI_Status stat;
                    //MPI::COMM_WORLD.Recv(&n_ex_ptcl_[tmp], 1, PS::GetDataType(*n_ex_ptcl_), i, TAG);
                    //MPI::COMM_WORLD.Recv(&n_ex_nei_[tmp],  1, PS::GetDataType(*n_ex_nei_),  i, TAG+1);
                    MPI_Recv(&n_ex_ptcl_[tmp], 1, PS::GetDataType(*n_ex_ptcl_), i, TAG,   MPI_COMM_WORLD, &stat);
                    MPI_Recv(&n_ex_nei_[tmp],  1, PS::GetDataType(*n_ex_nei_),  i, TAG+1, MPI_COMM_WORLD, &stat);
                    tmp ++;
                }
            }
        }
    }
}

template <class Tpsys>
void inputExParticle(Tpsys & pp,
                     FPGrav * & ex_ptcl_loc_,
                     PS::S32 * & ex_nei_loc_,
                     std::vector<std::vector<PS::S32> > & n_list,
                     const PS::S32 n_ex_ptcl_tot,
                     const PS::S32 n_ex_nei_tot)
{
    const PS::S32 n_loc = pp.getNumberOfParticleLocal();
    PS::S32 tmp = 0;
    PS::S32 tmp2 = 0;
    for(PS::S32 i=0; i<n_loc; ++i){
        if( pp[i].inDomain == false ){
            assert( pp[i].neighbor != 0 );
            for(PS::S32 j=0; j<pp[i].neighbor; ++j){
                ex_nei_loc_[tmp2] = n_list[i].at(j);
                tmp2 ++;
            }
            ex_ptcl_loc_[tmp] = pp[i];
            tmp ++;
        }
    }
    assert( tmp == n_ex_ptcl_tot );
    assert( tmp2 == n_ex_nei_tot );
}

template <class Tpsys>
void SendRecvOutOfDomainParticle(Tpsys & pp,
                                 PS::S32 * ex_rank,
                                 PS::S32 recv_rank,
                                 std::vector<std::vector<PS::S32> > n_list,
                                 PS::S32 ex_loc,
                                 PS::S32 ex_nei_loc,
                                 FPGrav * & ex_ptcl_,  
                                 PS::S32 * & ex_nei_,
                                 FPGrav * & ex_ptcl_loc_, 
                                 PS::S32 * & ex_nei_loc_,
                                 PS::S32 * & n_ex_ptcl_,
                                 PS::S32 * & n_ex_nei_)
{
    PS::S32 TAG = 2;
    if ( ex_loc ) {
        if ( PS::Comm::getRank() != recv_rank && ex_loc ){
            inputExParticle(pp, ex_ptcl_loc_, ex_nei_loc_, n_list, ex_loc, ex_nei_loc);
            //MPI::COMM_WORLD.Send(ex_ptcl_loc_, ex_loc,    PS::GetDataType(*ex_ptcl_loc_), recv_rank, TAG);
            //MPI::COMM_WORLD.Send(ex_nei_loc_,  ex_nei_loc,PS::GetDataType(*ex_nei_loc_),  recv_rank, TAG+1);
            MPI_Send(ex_ptcl_loc_, ex_loc,    PS::GetDataType(*ex_ptcl_loc_), recv_rank, TAG,   MPI_COMM_WORLD);
            MPI_Send(ex_nei_loc_,  ex_nei_loc,PS::GetDataType(*ex_nei_loc_),  recv_rank, TAG+1, MPI_COMM_WORLD);
        } else {
            PS::S32 nei_tmp = 0;
            PS::S32 nei_tmp2 = 0;
            PS::S32 tmp = 0;
            for(PS::S32 i=0; i<PS::Comm::getNumberOfProc(); ++i){
                if ( ex_rank[i] ){
                    if ( i == PS::Comm::getRank() ) continue;
                    MPI_Status stat;
                    //MPI::COMM_WORLD.Recv(&ex_ptcl_[nei_tmp], n_ex_ptcl_[tmp],PS::GetDataType(*ex_ptcl_), i, TAG);
                    //MPI::COMM_WORLD.Recv(&ex_nei_[nei_tmp2], n_ex_nei_[tmp], PS::GetDataType(*ex_nei_),  i, TAG+1);
                    MPI_Recv(&ex_ptcl_[nei_tmp], n_ex_ptcl_[tmp],PS::GetDataType(*ex_ptcl_), i, TAG,   MPI_COMM_WORLD, &stat);
                    MPI_Recv(&ex_nei_[nei_tmp2], n_ex_nei_[tmp], PS::GetDataType(*ex_nei_),  i, TAG+1, MPI_COMM_WORLD, &stat);
                    nei_tmp += n_ex_ptcl_[tmp];
                    nei_tmp2 += n_ex_nei_[tmp];
                    tmp ++;
                }
            }
        }
    }
}

template <class Tpsys>
void outputExParticle(Tpsys & pp,
                      FPGrav * & ex_ptcl_loc_,
                      const PS::S32 ex_loc)
{
    const PS::S32 n_loc = pp.getNumberOfParticleLocal();
    PS::S32  tmp = 0;
    for(PS::S32 i=0; i<n_loc; ++i){
        if( pp[i].inDomain == false ){
            assert( pp[i].neighbor != 0 );
            assert(ex_ptcl_loc_[tmp].id_local == pp[i].id_local);
            assert(ex_ptcl_loc_[tmp].myrank   == pp[i].myrank);
            pp[i] = ex_ptcl_loc_[tmp]; 
            tmp ++;
        }
    }
    assert(tmp == ex_loc);
}

template <class Tpsys>
void ReturnOutOfDomainParticle(Tpsys & pp,
                               PS::S32 * ex_rank,
                               PS::S32 recv_rank,
                               PS::S32 ex_loc,
                               FPGrav * ex_ptcl_,
                               FPGrav * ex_ptcl_loc_,
                               PS::S32 * n_ex_ptcl_)
{
    PS::S32 TAG = 4;
    if ( ex_loc ) {
        if ( PS::Comm::getRank() == recv_rank ){
            PS::S32 nei_tmp = 0;
            PS::S32 tmp = 0;
            for(PS::S32 i=0; i<PS::Comm::getNumberOfProc(); ++i){
                if ( ex_rank[i] ){
                    if ( i == PS::Comm::getRank() ) continue;
                    //MPI::COMM_WORLD.Send(ex_ptcl_+nei_tmp, n_ex_ptcl_[tmp], PS::GetDataType(*ex_ptcl_), i, TAG);
                    MPI_Send(ex_ptcl_+nei_tmp, n_ex_ptcl_[tmp], PS::GetDataType(*ex_ptcl_), i, TAG, MPI_COMM_WORLD);
                    nei_tmp += n_ex_ptcl_[tmp];
                    tmp ++;
                }
            }

        } else {
            MPI_Status stat;
            //MPI::COMM_WORLD.Recv(ex_ptcl_loc_, ex_loc, PS::GetDataType(*ex_ptcl_loc_), recv_rank, TAG);
            MPI_Recv(ex_ptcl_loc_, ex_loc, PS::GetDataType(*ex_ptcl_loc_), recv_rank, TAG, MPI_COMM_WORLD, &stat);
            outputExParticle(pp, ex_ptcl_loc_, ex_loc);
        }
    }
}
#endif

template <class Tpsys>
void MergeParticle(Tpsys & pp,
                   PS::S32 n_col,
                   PS::F64 & edisp)
{
    const PS::S32 n_loc = pp.getNumberOfParticleLocal();
    PS::S32 n_remove = 0;
    PS::S32 * remove = new PS::S32[n_col];
    PS::F64 edisp_loc = 0.;

#pragma omp parallel for reduction (-:edisp_loc)
    for ( PS::S32 i=0; i<n_loc; i++ ){
        if ( pp[i].isMerged ) {
            for ( PS::S32 j=0; j<n_loc; j++ ){              
                if ( pp[j].id == pp[i].id && i != j ){
                    PS::F64 mi = pp[i].mass;
                    PS::F64 mj = pp[j].mass;
                    PS::F64vec vrel = pp[j].vel - pp[i].vel;
                    pp[i].mass += mj;
                    pp[i].vel = ( mi*pp[i].vel + mj*pp[j].vel )/(mi+mj);
                    //pp[i].acc = ( mi*pp[i].acc + mj*pp[j].acc )/(mi+mj);
#ifdef GAS_DRAG
                    pp[i].acc_gd = ( mi*pp[i].acc_gd + mj*pp[j].acc_gd )/(mi+mj);
#endif
                    pp[i].phi   = ( mi*pp[i].phi   + mj*pp[j].phi   )/(mi+mj);
                    pp[i].phi_d = ( mi*pp[i].phi_d + mj*pp[j].phi_d )/(mi+mj);
                    
                    edisp_loc -= 0.5 * mi*mj/(mi+mj) * vrel*vrel;
#pragma omp critical
                    {
                        remove[n_remove] = j;
                        n_remove ++;
                    }
                    assert ( pp[i].pos == pp[j].pos );
                    assert ( pp[j].isDead );
                }
            }
            pp[i].isMerged = false;   
        }
    }
    PS::Comm::barrier();
    edisp += PS::Comm::getSum(edisp_loc);
    
    if ( n_remove ){
        pp.removeParticle(remove, n_remove);
    }
    delete [] remove;
}

template <class Tpsys>
void removeOutOfBoundaryParticle(Tpsys & pp,
                                 PS::F64 & edisp,
                                 const PS::F64 r_max,
                                 const PS::F64 r_min)
{
    PS::F64 rmax2 = r_max*r_max;
    PS::F64 rmin2 = r_min*r_min;
    PS::F64 edisp_loc = 0.;
    const PS::S32 n_loc = pp.getNumberOfParticleLocal();

    PS::S32 i_remove = -1;
    
#pragma omp parallel for
    for ( PS::S32 i=0; i<n_loc; i++ ){
        PS::F64vec posi = pp[i].pos;
        PS::F64    pos2 = posi*posi;
        if ( pos2 > rmax2 || pos2 < rmin2 ){
#pragma omp critical
            {
                i_remove = i;
            }
        }
    }

    if ( i_remove > -1 ){
        PS::F64    massi = pp[i_remove].mass;
        PS::F64vec veli = pp[i_remove].vel;
        edisp_loc -= 0.5*massi* veli*veli;
        edisp_loc -= massi * pp[i_remove].phi_s;
        edisp_loc -= massi * pp[i_remove].phi_d;
        edisp_loc -= massi * pp[i_remove].phi;

        pp.removeParticle(&i_remove, 1);
    }
    edisp += PS::Comm::getSum(edisp_loc);
}

template <class Tpsys>
void correctEnergyForGas(Tpsys & pp,
                         PS::F64 & edisp_gd,
                         bool second)
{// energy correction for gas drag
    PS::F64 edisp_gd_loc = 0.;
    PS::F64 coef = 0.25; if (second) coef *= -1.;
    const PS::S32 n_loc = pp.getNumberOfParticleLocal();
    
#pragma omp parallel for reduction(+:edisp_gd_loc)
    for(PS::S32 i=0; i<n_loc; i++){
        edisp_gd_loc += pp[i].mass * pp[i].acc_gd
            * (pp[i].vel + coef * pp[i].acc_gd * FPGrav::dt_tree);
    }
    edisp_gd += 0.5 * FPGrav::dt_tree * PS::Comm::getSum(edisp_gd_loc);
}
