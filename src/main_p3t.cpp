#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <stdlib.h>
#include <unistd.h>
#include <sys/stat.h>
#include <map>
#include <dirent.h>
#include <ctime>
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
#include <mpi.h>
#endif

#include <particle_simulator.hpp>

#define PRC(x) std::cerr << #x << " = " << x << ", "
#define PRL(x) std::cerr << #x << " = " << x << "\n"

#include "vector_x86.hpp"
#include "kepler.h"
#include "energy.h"
#include "particle.h"
#include "cutfunc.h"
#include "disk.h"
#include "gravity.h"
#include "collisionA.h"
#include "collisionB.h"
#include "hermite.h"
#include "hard.h"
#include "read.h"
#include "func.h"


int main(int argc, char *argv[])
{
    PS::Initialize(argc, argv);
    time_t wtime_start_program = time(NULL);
    
    ////////////////////////
    /*   Set Parameters   */
    ////////////////////////
    // Set Default Parameters
    char param_file[256] = "parameter.dat";
    
    char init_file[256] = "INIT_3000.dat";
    char output_dir[256] = "OUTPUT";
    bool existsHeader = false;
    bool isRestart    = false;

    bool makeInit = false;

    PS::F64 coef_ema = 0.3;
    PS::S32 nx = (int)sqrt(PS::Comm::getNumberOfProc());
    while ( PS::Comm::getNumberOfProc()%nx != 0 ) nx++;
    PS::S32 ny = PS::Comm::getNumberOfProc()/nx;
    
    PS::F64 theta         = 0.5;
    PS::S32 n_leaf_limit  = 8;
    PS::S32 n_group_limit = 256;    
    
    PS::F64 t_end   = pow2(-2);
    PS::F64 dt_snap = pow2(-5);

    PS::F64 r_max = 40.;
    PS::F64 r_min = 0.1;

    PS::S32 seed = 1;
    PS::F64 wtime_max = 0.;
    
    // Read Parameter File
    opterr = 0;
    PS::S32 opt;
    char param_file_opt[256];
    char init_file_opt[256];
    PS::S32 seed_opt;
    bool opt_p = false, opt_r = false, opt_i = false, opt_s = false;
    while ((opt = getopt(argc, argv, "prise:")) != -1) {
        switch (opt) {
        case 'p':
            sprintf(param_file_opt,"%s",optarg);
            opt_p = true;
            break;
        case 'r':
            opt_r = true;
            break;
        case 'i':
            sprintf(init_file_opt,"%s",optarg);
            opt_i = true;
            break;
        case 's':
            seed_opt = std::atoi(optarg);
            opt_s = true;
            break;
        case 'e':
            wtime_max = std::atof(optarg)*60.*60.;
            break;
        default:
            std::cout << "Usage: "
                      <<  argv[0]
                      << " [-p argment] [-r] [-i argment] [-s argment] [-e argment] arg1 ..." << std::endl;
            break;
        }
    }
    if ( readParameter(param_file, init_file, existsHeader, output_dir, isRestart,  makeInit,
                       coef_ema, nx, ny,
                       theta, n_leaf_limit, n_group_limit,
                       t_end, dt_snap, r_max, r_min, seed) ){
        PS::Abort();
        return 0;
    }
    if (opt_p) sprintf(param_file,"%s",param_file_opt);
    if (opt_r) isRestart = true;
    if (opt_i) sprintf(init_file,"%s",init_file_opt);
    if (opt_s) seed = seed_opt;
    
    char dir_name[256];
    sprintf(dir_name,"./%s",output_dir);
    if ( makeOutputDirectory(dir_name) ){
        PS::Abort();
        return 0;
    }
    if ( isRestart ) {
        if ( getLastSnap(dir_name, init_file) ) {
            PS::Abort();
            return 0;
        }
        existsHeader = true;
        makeInit = false;
    }

    srand48(seed);

    PS::F64 time_sys = 0.;
    PS::S32 n_tot  = 0;
    PS::F64 de_max = 0.;
    PS::S32 istep = 0;
    PS::S32 isnap = 0;
    PS::S32 n_col_tot  = 0;
    PS::S32 n_frag_tot = 0;
    
    ////////////////////////////////
    /*   Create Particle System   */
    ////////////////////////////////
    // Soft System
    PS::ParticleSystem<FPGrav> system_grav;
    system_grav.initialize();
    PS::S32 n_loc = 0;
    Energy e_init, e_now;
#ifdef OUTPUT_DETAIL
    PS::F64 ekin_before = 0., ekin_after = 0., edisp_d = 0.;
#endif
    PS::S32 id_next = 0;
    std::vector<std::vector<PS::S32> > n_list;
    n_list.clear();

    if ( makeInit && !isRestart){
        //Make Initial condition
        SolidDisk::createInitialCondition(system_grav);
        istep = 0;
        isnap = 0;
        time_sys = 0.;
        sprintf(init_file, "NONE");
    } else {
        // Read Initial File
        if ( existsHeader ) {
            FileHeader header;
            PS::F64 dt_tree = FPGrav::dt_tree;
            system_grav.readParticleAscii(init_file, header);
            if ( PS::Comm::getRank() == 0 ){
                istep = (PS::S32)round(header.time/dt_tree);
                isnap = (PS::S32)round(header.time/dt_snap);
                e_init = header.e_init;
                e_now = header.e_now;
            }
            PS::Comm::barrier();
            PS::Comm::broadcast(&istep, 1);
            PS::Comm::broadcast(&isnap, 1);
            PS::Comm::broadcast(&e_init, 1);
            time_sys = istep*dt_tree;
        } else {
            system_grav.readParticleAscii(init_file);
            istep = 0;
            isnap = 0;
            time_sys = 0.;
        }
    }
    n_loc = system_grav.getNumberOfParticleLocal();
    n_tot = system_grav.getNumberOfParticleGlobal();
    for ( PS::S32 i=0; i<n_loc; i++ ){
        if ( system_grav[i].id > id_next ) id_next = system_grav[i].id;
        system_grav[i].time = time_sys;
    }
    id_next = PS::Comm::getMaxValue(id_next);
    id_next ++;
    setCutoffRadii(system_grav);
        
    // Hard System
    HardSystem system_hard;
    system_hard.clear();
    PS::S32 nei_dist = 0;
    PS::S32 nei_tot_loc = 0;
    PS::S32 * ex_rank = new PS::S32[PS::Comm::getNumberOfProc()];
    PS::S32 * ex_rank_tmp = new PS::S32[PS::Comm::getNumberOfProc()];
    PS::S32 n_largestcluster = 0;

    ////////////////////
    /*   Set Domain   */
    ////////////////////
    PS::DomainInfo dinfo;
    dinfo.initialize(coef_ema);
    dinfo.setNumberOfDomainMultiDimension(nx,ny);
    dinfo.collectSampleParticle(system_grav);
    dinfo.decomposeDomain();
    system_grav.exchangeParticle(dinfo);
    n_loc = system_grav.getNumberOfParticleLocal();
    
    inputIDLocalAndMyrank(system_grav);
    
    /////////////////////
    /*   Create Tree   */
    /////////////////////
#ifdef USE_QUAD
    PS::TreeForForceLong<ForceGrav, EPGrav, EPGrav>::QuadrupoleWithScatterSearch tree_grav;
#else
    PS::TreeForForceLong<ForceGrav, EPGrav, EPGrav>::MonopoleWithScatterSearch tree_grav;
#endif
    tree_grav.initialize(n_tot, theta, n_leaf_limit, n_group_limit);
    tree_grav.calcForceAllAndWriteBack(CalcForceLongEP<EPGrav>,
#ifdef USE_QUAD
                                       CalcForceLongSP<PS::SPJQuadrupoleScatter>,
#else
                                       CalcForceLongSP<PS::SPJMonopoleScatter>,
#endif
                                       system_grav,
                                       dinfo);
    correctForceLongInitial(system_grav, tree_grav, n_list, nei_dist, nei_tot_loc, ex_rank);

#ifdef GAS_DRAG
    GasDisk gas_disk;
    gas_disk.calcGasDrag(system_grav, time_sys);
#endif
    
    //////////////////
    /*  File Open   */
    //////////////////
    std::ofstream fout_eng;
    std::ofstream fout_col;
    char sout_eng[256];
    char sout_col[256];
    sprintf(sout_eng, "%s/energy.dat", dir_name);
    sprintf(sout_col, "%s/collision.dat", dir_name);
    if ( time_sys == 0. ) {
        fout_eng.open(sout_eng, std::ios::out);
        fout_col.open(sout_col, std::ios::out);
    } else {
        fout_eng.open(sout_eng, std::ios::app);
        fout_col.open(sout_col, std::ios::app);
    }
    PS::Comm::barrier();
    showParameter(init_file, dir_name, makeInit,
                  time_sys,
                  coef_ema, nx, ny,
                  theta, n_leaf_limit, n_group_limit,
                  t_end, dt_snap, r_max, r_min, seed);

    ////////////////////////////////
    /*  Preparation Before Loop   */
    ////////////////////////////////
    e_now.calcEnergy(system_grav);
    if ( !existsHeader ) e_init = e_now;
    PS::F64 de =  e_now.calcEnergyError(e_init);
    
    PS::F64 wtime_init = 0.0;
    PS::F64 wtime_now = 0.0;
    PS::F64 wtime_soft = 0.0;
    PS::F64 wtime_hard = 0.0;
    PS::F64 wtime_soft_step = 0.0;
    PS::F64 wtime_hard_step = 0.0;
    PS::F64 wtime_start_soft = 0.0;
    PS::F64 wtime_end_soft = 0.0;
    PS::F64 wtime_start_hard = 0.0;
    PS::F64 wtime_end_hard = 0.0;  

    PS::Comm::barrier();
    wtime_init = wtime_now = PS::GetWtime();

    outputStep(system_grav, time_sys, e_init, e_now, de,
               n_col_tot, n_frag_tot, dir_name, isnap, fout_eng, 
               wtime_soft_step, wtime_hard_step, n_largestcluster,
               (time_sys==0.) );
    istep ++;
    isnap ++;
    
    //////////////////////
    ///   Loop Start
    while(time_sys < t_end){
        
        PS::S32 n_col = 0;
        PS::S32 n_frag = 0;
        std::map<PS::S32, PS::S32> id2id_loc;
        id2id_loc.clear();
#ifdef GAS_DRAG
        PS::F64 edisp_gd = 0.;
#endif
        
        n_loc = system_grav.getNumberOfParticleLocal();

        PS::Comm::barrier();
        wtime_start_soft = PS::GetWtime();

        ////////////////////
        ///   Soft Part
        
        ///////////////////////////
        /*   1st Velocity kick   */
        ///////////////////////////
#ifdef GAS_DRAG
        correctEnergyForGas(system_grav, edisp_gd, false);
#endif
        velKick(system_grav);
#ifdef OUTPUT_DETAIL
        calcKineticEnergy(system_grav, ekin_before);
#endif
        
        ///   Soft Part
        ////////////////////
        

        PS::Comm::barrier();
        wtime_end_soft = wtime_start_hard = PS::GetWtime();
        wtime_soft += wtime_soft_step = wtime_end_soft - wtime_start_soft;
        
        
        /////////////////////
        ///   Hard Part

        ////////////////////////////
        /*   Create Hard System   */
        ////////////////////////////
        PS::S32 n_ex_rank = 0;
        PS::S32 recv_rank = PS::Comm::getRank();
        PS::S32 ex_tot = 0;
        PS::S32 ex_loc = 0;
        PS::S32 ex_nei_loc = 0;

        createNeighborCluster(system_grav, id2id_loc, n_list);

        decideExchangeRank(ex_rank, ex_rank_tmp, n_ex_rank, recv_rank); 
        countSendParticle(system_grav, ex_tot, ex_loc, ex_nei_loc); 
        //std::cerr << "Rank:" << PS::Comm::getRank() << " RecvRank:" << recv_rank
        //          << " ex_loc = " << ex_loc << std::endl;

        FPGrav * ex_ptcl_ = nullptr;
        PS::S32 * ex_nei_ = nullptr;
        FPGrav * ex_ptcl_loc_ = nullptr;
        PS::S32 * ex_nei_loc_ = nullptr;
        PS::S32 * n_ex_ptcl_ = nullptr;
        PS::S32 * n_ex_nei_  = nullptr;
        PS::S32 * ex_nei_adr_ = nullptr;
        PS::S32 n_ex_ptcl_tot = 0;
        PS::S32 n_ex_nei_tot  = 0;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        if( ex_tot != 0 ){            
            /////////////////////
            /*   Gather Hard   */
            /////////////////////
            if ( PS::Comm::getRank() == recv_rank && ex_loc ){
                n_ex_ptcl_ = new PS::S32[n_ex_rank -1];
                n_ex_nei_  = new PS::S32[n_ex_rank -1];
            }
            //PS::Comm::barrier();

            // Send & Receive Particles Number
            SendRecvOutOfDomainParticleNumber(ex_rank, recv_rank, ex_loc, ex_nei_loc,
                                              n_ex_ptcl_, n_ex_nei_);

            if ( PS::Comm::getRank() == recv_rank && ex_loc ){
                for ( PS::S32 i=0; i<n_ex_rank-1; i++ ){
                    n_ex_ptcl_tot += n_ex_ptcl_[i];
                    n_ex_nei_tot  += n_ex_nei_ [i];
                }
                ex_ptcl_ = new FPGrav[n_ex_ptcl_tot];
                ex_nei_  = new PS::S32[n_ex_nei_tot];
                ex_nei_adr_ = new PS::S32[n_ex_ptcl_tot];
            } else if ( ex_loc ){
                ex_ptcl_loc_ = new FPGrav[ex_loc];
                ex_nei_loc_  = new PS::S32[ex_nei_loc];
            }
            //PS::Comm::barrier();

            // Send & Receive Particles
            SendRecvOutOfDomainParticle(system_grav, ex_rank, recv_rank, n_list,
                                        ex_loc, ex_nei_loc,
                                        ex_ptcl_, ex_nei_, ex_ptcl_loc_, ex_nei_loc_,
                                        n_ex_ptcl_, n_ex_nei_);

            if ( PS::Comm::getRank() == recv_rank && ex_loc ){
                ex_nei_adr_[0] = 0;
                for ( PS::S32 i=1; i<n_ex_ptcl_tot; i++ ){
                    ex_nei_adr_[i] = ex_nei_adr_[i-1] + ex_ptcl_[i-1].neighbor;
                }
                assert( ex_nei_adr_[n_ex_ptcl_tot-1]+ex_ptcl_[n_ex_ptcl_tot-1].neighbor == n_ex_nei_tot);
                
                createNeighborCluster_OutOfDomain(system_grav, n_loc, ex_ptcl_, n_ex_ptcl_tot,
                                                  id2id_loc, n_list, ex_nei_, ex_nei_adr_, //n_ex_nei_tot,
                                                  recv_rank);
            }
        }
#endif

        ////////////////////////
        /*   Time Integrate   */
        ////////////////////////
        PS::S32 n_in = 0;
        PS::S32 n_out = 0;
        system_hard.clear();
        n_in = system_hard.makeList(system_grav, n_loc, ex_ptcl_, n_ex_ptcl_tot);
        n_out = system_hard.timeIntegrate(system_grav, ex_ptcl_, n_list, ex_nei_, ex_nei_adr_, istep);
        //if( ex_tot != 0 )
        //system_hard.showParticleID();
        assert( n_in == n_out );      
        
        //PS::Comm::barrier();
        n_col  = system_hard.getNumberOfCollisionGlobal();
        //n_frag = system_hard.getNumberOfFragmentGlobal();
        if ( n_col ) n_frag = system_hard.addFragment2ParticleSystem(system_grav, id_next, fout_col);
        
        n_col_tot   += n_col;
        n_frag_tot  += n_frag;
        e_now.edisp += system_hard.getEnergyDissipationGlobal();
#ifdef OUTPUT_DETAIL
        edisp_d      = system_hard.getHardEnergyDissipationGlobal();
#endif
        n_largestcluster = system_hard.getNumberOfParticleInLargestClusterGlobal();

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        if( ex_tot != 0 ){
            //////////////////////
            /*   Scatter Hard   */
            //////////////////////
            ReturnOutOfDomainParticle(system_grav, ex_rank, recv_rank, ex_loc,
                                      ex_ptcl_, ex_ptcl_loc_, n_ex_ptcl_);
            //PS::Comm::barrier();

            if ( PS::Comm::getRank() == recv_rank && ex_loc ){
                delete [] ex_ptcl_;
                delete [] ex_nei_;
                delete [] n_ex_ptcl_;
                delete [] n_ex_nei_;
                delete [] ex_nei_adr_;
            } else if ( ex_loc ) {
                delete [] ex_ptcl_loc_;
                delete [] ex_nei_loc_;
            }           
        }
#endif
        ///   Hard Part
        ////////////////////
          
   
        PS::Comm::barrier();
        wtime_end_hard = wtime_start_soft = PS::GetWtime();
        wtime_hard += wtime_hard_step = wtime_end_hard - wtime_start_hard;
        

        ////////////////////
        ///   Soft Part

        ////////////////////////
        /*   Calculate Soft   */
        ////////////////////////
        system_grav.exchangeParticle(dinfo);
        inputIDLocalAndMyrank(system_grav);
        tree_grav.calcForceAllAndWriteBack(CalcForceLongEP<EPGrav>,
#ifdef USE_QUAD
                                           CalcForceLongSP<PS::SPJQuadrupoleScatter>,
#else
                                           CalcForceLongSP<PS::SPJMonopoleScatter>,
#endif
                                           system_grav,
                                           dinfo);
        correctForceLong(system_grav, tree_grav, n_list, nei_dist, nei_tot_loc, ex_rank);
        
#ifdef GAS_DRAG
        gas_disk.calcGasDrag(system_grav, time_sys+FPGrav::dt_tree);
#endif

        ///////////////////////////
        /*   2nd Velocity kick   */
        ///////////////////////////
#ifdef OUTPUT_DETAIL
        calcKineticEnergy(system_grav, ekin_after);
#endif
        velKick(system_grav);
#ifdef GAS_DRAG
        correctEnergyForGas(system_grav, edisp_gd, true);
        e_now.edisp += edisp_gd;
#endif
        
        //////////////////////////
        /*   Calculate Energy   */
        //////////////////////////
#ifdef OUTPUT_DETAIL
        Energy e_tmp = e_now;
        PS::F64 dekin_d = ekin_after - ekin_before;
#endif
        e_now.calcEnergy(system_grav);
#ifdef OUTPUT_DETAIL
        PS::F64 dephi_d_d = e_now.ephi_d - e_tmp.ephi_d;
        PS::F64 dephi_s_d = e_now.ephi_sun - e_tmp.ephi_sun;
        PS::F64 de_d = ( dekin_d + dephi_d_d + dephi_s_d - edisp_d ) / e_init.etot;
#endif

        ///////////////
        /*   Merge   */
        ///////////////
        if ( n_col ) MergeParticle(system_grav, n_col, e_now.edisp);

        // Reset Number Of Particles
        n_tot = system_grav.getNumberOfParticleGlobal();
        n_loc = system_grav.getNumberOfParticleLocal();
        
        ///////////////////////////
        /*   Re-Calculate Soft   */
        ///////////////////////////
        if ( n_col || istep % 128 == 0 ) {
            if(istep % 128 == 0) {
                dinfo.decomposeDomainAll(system_grav);
                system_grav.exchangeParticle(dinfo);

                // Remove Particle Out Of Boundary
                removeOutOfBoundaryParticle(system_grav, e_now.edisp, r_max, r_min);
                
                // Reset Number Of Particles
                n_tot = system_grav.getNumberOfParticleGlobal();
                n_loc = system_grav.getNumberOfParticleLocal();
            }
            inputIDLocalAndMyrank(system_grav);
            
            setCutoffRadii(system_grav);
            tree_grav.calcForceAllAndWriteBack(CalcForceLongEP<EPGrav>,
#ifdef USE_QUAD
                                               CalcForceLongSP<PS::SPJQuadrupoleScatter>,
#else
                                               CalcForceLongSP<PS::SPJMonopoleScatter>,
#endif
                                               system_grav,
                                               dinfo);
            correctForceLongInitial(system_grav, tree_grav, n_list, nei_dist, nei_tot_loc, ex_rank);
#ifdef GAS_DRAG
#pragma omp parallel for
            for(PS::S32 i=0; i<n_loc; i++) system_grav[i].acc += system_grav[i].acc_gd;
#endif
            e_now.calcEnergy(system_grav);
        }
   
        ///   Soft Part
        ////////////////////

        PS::Comm::barrier();
        wtime_now = wtime_end_soft = PS::GetWtime();
        wtime_soft += wtime_end_soft - wtime_start_soft;
        wtime_soft_step +=  wtime_end_soft - wtime_start_soft;
        
        ////////////////
        /*   Output   */
        ////////////////
        time_sys += FPGrav::dt_tree;
        istep ++;

        PS::F64 de =  e_now.calcEnergyError(e_init);
        PS::F64 de_tmp = sqrt(de*de);
        if( de_tmp > de_max ) de_max = de_tmp;
        if ( PS::Comm::getRank() == 0 ) {
            std::cout << std::fixed << std::setprecision(8)
                      << "Time: " << time_sys
                      << "  NumberOfParticle: " << n_tot 
                      << "  NumberOfCol: " << n_col_tot
                      << "  NumberOfFrag: " << n_frag_tot << std::endl
                      << std::scientific << std::setprecision(15)
                      << "                  "
#ifdef OUTPUT_DETAIL
                      << "HardEnergyError: " << de_d
                      << "  "
#endif
                      << "EnergyError: " << de
                      << "  MaxEnergyError: " << de_max << std::endl;
            
            std::cerr << std::scientific<<std::setprecision(15);
            //PRC(etot1); PRL(ekin);
            //PRC(ephi); PRC(ephi_s); PRL(ephi_d);
        }
        
        if( time_sys  == dt_snap*isnap ){
            outputStep(system_grav, time_sys, e_init, e_now, de,
                       n_col_tot, n_frag_tot, dir_name, isnap, fout_eng,
                       wtime_soft_step, wtime_hard_step, n_largestcluster);
            isnap ++;

            if ( wtime_max > 0. && wtime_max < difftime(time(NULL), wtime_start_program) ) break; 
        }
    }
       
    ///   Loop End
    ////////////////////
    
    delete[] ex_rank;
    delete[] ex_rank_tmp;
    fout_eng.close();
    fout_col.close();
    
    PS::Comm::barrier();

    wtime_now = PS::GetWtime();
    showTime(dir_name, wtime_start_program, wtime_init, wtime_now, wtime_soft, wtime_hard);
    
    PS::Finalize();
    return 0;
}
