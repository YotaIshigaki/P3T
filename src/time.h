#pragma once

class Wtime{
public:
    PS::F64 init;
    PS::F64 now;

    PS::F64 soft;
    PS::F64 hard;
    PS::F64 soft_step;
    PS::F64 hard_step;

    PS::F64 start_soft;
    PS::F64 end_soft;
    PS::F64 start_hard;
    PS::F64 end_hard;

#ifdef CALC_WTIME
    PS::F64 calc_soft_force;
    PS::F64 calc_hard_force;
    PS::F64 create_cluster;
    PS::F64 communication;
    PS::F64 output;
    PS::F64 calc_soft_force_step;
    PS::F64 calc_hard_force_step;
    PS::F64 create_cluster_step;
    PS::F64 communication_step;
    PS::F64 output_step;

    PS::F64 laptime;

    PS::F64 lap(PS::F64 time_now){
        PS::F64 time = time_now - laptime;
        laptime = time_now;
        return time;
    }
#endif

    Wtime(){
        init = now = soft = hard = 0.;
        start_soft = end_soft = start_hard = end_hard = 0.;
#ifdef CALC_WTIME
        calc_soft_force = output = calc_hard_force
            = create_cluster = communication = laptime = 0.;
#endif
    }

    void showTime(char * dir_name,
                  time_t wtime_start_program);
};

void Wtime::showTime(char * dir_name,
              time_t wtime_start_program)
{
    if( PS::Comm::getRank() == 0 ){
        char wtime_s[64], wtime_n[64];
        time_t wtime_now_program = time(NULL);
        strftime(wtime_s, sizeof(wtime_s), "%Y/%m/%d %a %H:%M:%S", localtime(&wtime_start_program));
        strftime(wtime_n, sizeof(wtime_n), "%Y/%m/%d %a %H:%M:%S", localtime(&wtime_now_program));      
        
        char sout_param[256];
        std::ofstream fout_param;
        sprintf(sout_param,"%s/param.dat", dir_name);
        fout_param.open(sout_param, std::ios::app);
        fout_param << std::endl
                   << "Start Wall Time: " << wtime_s << std::endl
                   << "End Wall Time:   " << wtime_n << std::endl
                   << std::scientific << std::setprecision(15)
                   << "Wall Time: " << now - init
                   << "\tSoft: " << soft
                   << "\tHard: " << hard << std::endl;
        
        std::cout << "Wall Time: " << now - init
                  << "\tSoft: " << soft
                  << "\tHard: " << hard << std::endl;

#ifdef CALC_WTIME
        fout_param << "\tCalc Soft Force: " << calc_soft_force
                   << "\tCalc Hard Force: " <<  calc_hard_force
                   << "\tCreate Cluster: " << create_cluster
                   << "\tCommunication: " << communication << std::endl
                   << "\tOutput: " << output << std::endl;

        std::cout << "\tCalc Soft Force: " << calc_soft_force
                  << "\tCalc Hard Force: " <<  calc_hard_force
                  << "\tCreate Cluster: " << create_cluster
                  << "\tCommunication: " << communication << std::endl
                  << "\tOutput: " << output << std::endl;
#endif
    }
}
