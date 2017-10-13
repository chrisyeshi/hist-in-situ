#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#include <stdio.h>
#include <stdint.h>
#include <mpi.h>

void c_generate_and_output_histogram_1d_(
        double *ptr_values,
        double *ptr_log_base,
        char* ptr_method,
        double *ptr_min, double *ptr_max,
        double *ptr_timestep,
        int32_t *ptr_ngridx, int32_t *ptr_ngridy, int32_t *ptr_ngridz,
        int32_t *ptr_nhistx, int32_t *ptr_nhisty, int32_t *ptr_nhistz,
        char *ptr_nbin,
        int32_t *ptr_yid, int32_t *ptr_rootid, int32_t *ptr_ypes,
        MPI_Fint *ptr_comm, int32_t *ptr_xzid, int32_t *ptr_iconfig);

void c_generate_and_output_histogram_2d_(
        double *ptr_values_x, double *ptr_values_y,
        double *ptr_log_base_x, double *ptr_log_base_y,
        char *ptr_method_x, char *ptr_method_y,
        double *ptr_min_x, double *ptr_min_y,
        double *ptr_max_x, double *ptr_max_y,
        double *ptr_timestep,
        int32_t *ptr_ngridx, int32_t *ptr_ngridy, int32_t *ptr_ngridz,
        int32_t *ptr_nhistx, int32_t *ptr_nhisty, int32_t *ptr_nhistz,
        char *ptr_nbinx, char *ptr_nbiny,
        int32_t *ptr_yid, int32_t *ptr_rootid, int32_t *ptr_ypes,
        MPI_Fint *ptr_comm, int32_t *ptr_xzid, int32_t *ptr_iconfig);

void c_generate_and_output_histogram_3d_(
        double *ptr_values_x, double *ptr_values_y, double *ptr_values_z,
        double *ptr_log_base_x, double *ptr_log_base_y, double *ptr_log_base_z,
        char *ptr_method_x, char *ptr_method_y, char* ptr_method_z,
        double *ptr_min_x, double *ptr_min_y, double *ptr_min_z,
        double *ptr_max_x, double *ptr_max_y, double *ptr_max_z,
        double *ptr_timestep,
        int32_t *ptr_ngridx, int32_t *ptr_ngridy, int32_t *ptr_ngridz,
        int32_t *ptr_nhistx, int32_t *ptr_nhisty, int32_t *ptr_nhistz,
        char *ptr_nbinx, char *ptr_nbiny, char *ptr_nbinz,
        int32_t *ptr_yid, int32_t *ptr_rootid, int32_t *ptr_ypes,
        MPI_Fint *ptr_comm, int32_t *ptr_xzid, int32_t *ptr_iconfig);

void c_generate_and_output_histogram_(
        double *ptr_values, double *ptr_log_bases, char* ptr_methods,
        char *ptr_nbins, double *ptr_mins, double *ptr_maxs,
        double *ptr_timestep, int32_t *ptr_ndims,
        int32_t *ptr_ngridx, int32_t *ptr_ngridy, int32_t *ptr_ngridz,
        int32_t *ptr_nhistx, int32_t *ptr_nhisty, int32_t *ptr_nhistz,
        int32_t *ptr_yid, int32_t *ptr_rootid,
        int32_t *ptr_ypes, MPI_Fint *ptr_comm, int32_t *ptr_xzid,
        int32_t *ptr_iconfig);

#endif // HISTOGRAM_H
