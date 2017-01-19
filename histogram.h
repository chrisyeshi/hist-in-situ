#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#include <stdio.h>
#include <stdint.h>
#include <mpi.h>

#define MAX_DIM 3

typedef struct {
    int32_t ndims;
    double* valuearrays[MAX_DIM];
    int32_t ngridx, ngridy, ngridz;
    int32_t xbeg, xend, ybeg, yend, zbeg, zend;
} SamplingRegion;

typedef struct {
    int32_t ndims;
    int32_t* buffer;
    int issparse;
    int nbins[MAX_DIM];
    int nnonemptybins;
    double percentinrange;
    double mins[MAX_DIM], maxs[MAX_DIM];
} Histogram;

void c_generate_and_output_histogram_1d_(
        double *ptr_values,
        double *ptr_log_base,
        char* ptr_method,
        double *ptr_min, double *ptr_max,
        double *ptr_timestep,
        int32_t *ptr_ngridx, int32_t *ptr_ngridy, int32_t *ptr_ngridz,
        int32_t *ptr_nhistx, int32_t *ptr_nhisty, int32_t *ptr_nhistz,
        int32_t *ptr_nbin,
        int32_t *ptr_myid, int32_t *ptr_rootid, MPI_Fint *ptr_comm,
        int32_t *ptr_iconfig);

void c_generate_and_output_histogram_2d_(
        double *ptr_values_x, double *ptr_values_y,
        double *ptr_log_base_x, double *ptr_log_base_y,
        char *ptr_method_x, char *ptr_method_y,
        double *ptr_min_x, double *ptr_min_y,
        double *ptr_max_x, double *ptr_max_y,
        double *ptr_timestep,
        int32_t *ptr_ngridx, int32_t *ptr_ngridy, int32_t *ptr_ngridz,
        int32_t *ptr_nhistx, int32_t *ptr_nhisty, int32_t *ptr_nhistz,
        int32_t *ptr_nbinx, int32_t *ptr_nbiny,
        int32_t *ptr_myid, int32_t *ptr_rootid, MPI_Fint *ptr_comm,
        int32_t *ptr_iconfig);

void c_generate_and_output_histogram_3d_(
        double *ptr_values_x, double *ptr_values_y, double *ptr_values_z,
        double *ptr_log_base_x, double *ptr_log_base_y, double *ptr_log_base_z,
        char *ptr_method_x, char *ptr_method_y, char* ptr_method_z,
        double *ptr_min_x, double *ptr_min_y, double *ptr_min_z,
        double *ptr_max_x, double *ptr_max_y, double *ptr_max_z,
        double *ptr_timestep,
        int32_t *ptr_ngridx, int32_t *ptr_ngridy, int32_t *ptr_ngridz,
        int32_t *ptr_nhistx, int32_t *ptr_nhisty, int32_t *ptr_nhistz,
        int32_t *ptr_nbinx, int32_t *ptr_nbiny, int32_t *ptr_nbinz,
        int32_t *ptr_myid, int32_t *ptr_rootid, MPI_Fint *ptr_comm,
        int32_t *ptr_iconfig);

int32_t total_number_of_bins(int32_t ndims, int32_t nbins[]);

int32_t get_number_of_nonempty_bins(int32_t* frequencies, int32_t nbins);

Histogram frequencies_to_histogram_dense(
        int32_t ndims, int32_t* frequencies, int32_t nbins[],
        int32_t nnonemptybins, double mins[], double maxs[],
        double percentinrange);

Histogram frequencies_to_histogram_sparse(
        int32_t ndims, int32_t* frequencies, int32_t nbins[],
        int32_t nnonemptybins, double mins[], double maxs[],
        double percentinrange);

Histogram frequencies_to_histogram(
        int32_t ndims, int32_t* frequencies, int32_t nbins[],
        double mins[], double maxs[], double percentinrange);

void adjust_range_by_methods(int32_t ndims, SamplingRegion samplingregion,
        char* methods[], double mins[], double maxs[], double adjmins[],
        double adjmaxs[]);

Histogram generate_histogram_from_sampling_region(
        int32_t ndims, SamplingRegion samplingregion, int32_t nbins[],
        double mins[], double maxs[]);

void deallocate_histogram(Histogram hist);

void print_histogram_meta(
        int32_t ngridx, int32_t ngridy, int32_t ngridz,
        int32_t nhistx, int32_t nhisty, int32_t nhistz,
        int32_t ndims, int32_t nbins[], double log_base[]);

void print_histogram(Histogram hist);

void print_domain_histograms(
        int32_t ngridx, int32_t ngridy, int32_t ngridz,
        int32_t nhistx, int32_t nhisty, int32_t nhistz,
        int32_t ndims, int32_t nbins[], double log_base[],
        Histogram hists[], int32_t nhists);

void write_histogram_meta(
        FILE* outfile,
        int32_t ngridx, int32_t ngridy, int32_t ngridz,
        int32_t nhistx, int32_t nhisty, int32_t nhistz,
        int32_t ndims, int32_t nbins[], double log_base[]);

void write_histogram(FILE* outfile, Histogram hist);

void write_domain_histograms(
        int32_t ngridx, int32_t ngridy, int32_t ngridz,
        int32_t nhistx, int32_t nhisty, int32_t nhistz,
        int32_t ndims, int32_t nbins[], double log_base[],
        Histogram hists[], int32_t nhists,
        double timestep, int32_t myid, int32_t iconfig);

SamplingRegion construct_sampling_region(
        double* valuearrays[], int32_t ndims,
        int32_t ngridx, int32_t ngridy, int32_t ngridz,
        int32_t nhistx, int32_t nhisty, int32_t nhistz,
        int32_t ihistx, int32_t ihisty, int32_t ihistz);

int32_t ids_to_flat(int32_t ndims, int32_t dims[], int32_t ids[]);

int32_t values_to_bin_index(int32_t ndims, double values[],
        double mins[], double maxs[], int32_t nbins[]);

void generate_and_output_histogram(
        double* values[], double log_bases[], char* methods[],
        double mins[], double maxs[],
        double timestep,
        int32_t ndims,
        int32_t ngridx, int32_t ngridy, int32_t ngridz,
        int32_t nhistx, int32_t nhisty, int32_t nhistz,
        int32_t nbins[], int32_t myid, int32_t iconfig);

#endif // HISTOGRAM_H
