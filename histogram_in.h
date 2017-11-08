#ifndef HISTOGRAM_IN_H
#define HISTOGRAM_IN_H

#define MAX_DIM 3
#define MAX_NBINS 100

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

void compute_nbins(int32_t ndims, SamplingRegion samplingregion, double mins[],
        double maxs[], char* nbinstrs[], int32_t nbins[]);

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

void serialize_histogram_meta(
        int32_t ngridx, int32_t ngridy, int32_t ngridz,
        int32_t nhistx, int32_t nhisty, int32_t nhistz,
        int32_t ndims, double log_base[], char** buffer, int32_t* nbyte);

void deserialize_histogram_meta(
        char* buffer, int32_t nbyte,
        int32_t* ngridx, int32_t* ngridy, int32_t* ngridz,
        int32_t* nhistx, int32_t* nhisty, int32_t* nhistz,
        int32_t* ndims, double log_base[]);

void write_histogram(FILE* outfile, Histogram hist);

int hist_buffer_byte_count(int32_t issparse, int32_t ndims, int32_t nbins[],
        int32_t nnonemptybins);

int hist_byte_count(int32_t issparse, int32_t ndims, int32_t nbins[],
        int32_t nnonemptybins);

int buffer_to_hist_byte(char* buffer, int ndims);

int buffer_to_hist_count(int nBytes, char* buffer, int ndims);

void serialize_histogram(Histogram hist, char** buffer, int32_t *nbyte);

void serialize_histograms(
        int nHists, Histogram hists[], int* nBytes, char** buffer);

void deserialize_histogram(
        char* buffer, int ndims, int* histBufferByteCount, Histogram* hist);

void deserialize_histograms(
        int nBytes, char* buffer, int ndims, int* nHists, Histogram** hists);

int is_same_histogram(Histogram a, Histogram b);

void write_domain_histograms(
        int32_t ngridx, int32_t ngridy, int32_t ngridz,
        int32_t nhistx, int32_t nhisty, int32_t nhistz,
        int32_t ndims, int32_t nbins[], double log_base[],
        Histogram hists[], int32_t nhists,
        double timestep, int32_t yid, int32_t rootid, MPI_Comm comm,
        int32_t xzid, int32_t iconfig);

void serialize_domain_histograms(
        int32_t ngridx, int32_t ngridy, int32_t ngridz,
        int32_t nhistx, int32_t nhisty, int32_t nhistz,
        int32_t ndims, double log_base[],
        Histogram hists[], int32_t nhists, char** buffer, int32_t* nbyte);

void write_histograms_grouped_by_mpi_comm(
        int32_t ngridx, int32_t ngridy, int32_t ngridz,
        int32_t nhistx, int32_t nhisty, int32_t nhistz,
        int32_t ndims, double log_base[], Histogram hists[], int32_t nhists,
        double timestep, int32_t yid, int32_t rootid, int32_t ypes,
        MPI_Comm comm, int32_t xzid, int32_t iconfig, const char* midname);

void write_ycolumn_histograms(
        int32_t ngridx, int32_t ngridy, int32_t ngridz,
        int32_t nhistx, int32_t nhisty, int32_t nhistz,
        int32_t ndims, double log_base[], Histogram hists[], int32_t nhists,
        double timestep, int32_t yid, int32_t rootid, int32_t ypes,
        MPI_Comm comm, int32_t xzid, int32_t iconfig);

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
        char* nbins[], int32_t myid, int32_t rootid, int32_t ypes,
        MPI_Comm comm, int32_t xzid, int32_t iconfig);

//
// mpi helper functions
//

void gather_buffers(int nBytesPerRank, char* bufferPerRank,
        int** nBytesPerRanks, int* nBytes, char** buffer, int root,
        MPI_Comm comm);

#endif // HISTOGRAM_IN_H