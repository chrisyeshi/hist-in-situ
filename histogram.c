#include "histogram.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <float.h>
#include <assert.h>
#include <math.h>
#include <mpi.h>

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

/// TODO: change histogram dimension to a, b, c

#define PRINT_HISTOGRAM

int32_t total_number_of_bins(int32_t ndims, int32_t nbins[]) {
    int32_t i, total = 1;
    for (i = 0; i < ndims; ++i)
        total *= nbins[i];
    return total;
}

int32_t get_number_of_nonempty_bins(int32_t* frequencies, int32_t nbins) {
    int32_t nnonemptybins = 0, ibin;
    assert(nbins >= 0);
    for (ibin = 0; ibin < nbins; ++ibin) {
        if (frequencies[ibin] > 0)
            ++nnonemptybins;
    }
    return nnonemptybins;
}

Histogram frequencies_to_histogram_dense(
        int32_t ndims, int32_t* frequencies, int32_t nbins[], int32_t nnonemptybins,
        double mins[], double maxs[], double percentinrange) {
    int32_t i, totalbin;
    Histogram hist;
    hist.ndims = ndims;
    hist.issparse = 0;
    hist.percentinrange = percentinrange;
    hist.nnonemptybins = nnonemptybins;
    for (i = 0; i < ndims; ++i) {
        hist.nbins[i] = nbins[i];
        hist.mins[i] = mins[i];
        hist.maxs[i] = maxs[i];
    }
    totalbin = total_number_of_bins(ndims, nbins);
    hist.buffer = malloc(totalbin * sizeof(int32_t));
    memcpy(hist.buffer, frequencies, totalbin * sizeof(int32_t));
    return hist;
}

Histogram frequencies_to_histogram_sparse(
        int32_t ndims, int32_t* frequencies, int32_t nbins[], int32_t nnonemptybins,
        double mins[], double maxs[], double percentinrange) {
    int32_t i, ibin, sparseindex, totalbin;
    Histogram hist;
    hist.ndims = ndims;
    hist.issparse = 1;
    hist.percentinrange = percentinrange;
    hist.nnonemptybins = nnonemptybins;
    for (i = 0; i < ndims; ++i) {
        hist.nbins[i] = nbins[i];
        hist.mins[i] = mins[i];
        hist.maxs[i] = maxs[i];
    }
    hist.buffer = malloc(2 * nnonemptybins * sizeof(int32_t));
    sparseindex = 0;
    totalbin = total_number_of_bins(ndims, nbins);
    for (ibin = 0; ibin < totalbin; ++ibin) {
        if (frequencies[ibin] > 0) {
            hist.buffer[sparseindex++] = ibin;
            hist.buffer[sparseindex++] = frequencies[ibin];
        }
    }
    return hist;
}

Histogram frequencies_to_histogram(
        int32_t ndims, int32_t* frequencies, int32_t nbins[],
        double mins[], double maxs[], double percentinrange) {
    int32_t nnonemptybins, i, ibin, nbin;
    nbin = total_number_of_bins(ndims, nbins);
    nnonemptybins = get_number_of_nonempty_bins(frequencies, nbin);
    if (nnonemptybins > 0.5 * nbin) {
        return frequencies_to_histogram_dense(
                ndims, frequencies, nbins, nnonemptybins, mins, maxs, percentinrange);
    }
    return frequencies_to_histogram_sparse(
            ndims, frequencies, nbins, nnonemptybins, mins, maxs, percentinrange);
}

void deallocate_histogram(Histogram hist) {
    free(hist.buffer);
}

void print_histogram_meta(
        int32_t ngridx, int32_t ngridy, int32_t ngridz,
        int32_t nhistx, int32_t nhisty, int32_t nhistz,
        int32_t ndims, int32_t nbins[], double log_base[]) {
    int32_t i;
    printf("ndims:%d    ", ndims);
    printf("nvoxels:%d %d %d    ", ngridx, ngridy, ngridz);
    printf("nhists:%d %d %d    ", nhistx, nhisty, nhistz);
    printf("nbins:");
    for (i = 0; i < ndims; ++i)
        printf("%d ", nbins[i]);
    printf("   ");
    printf("logbase:");
    for (i = 0; i < ndims; ++i)
        printf("%lf ", log_base[i]);
    printf("   \n");
}

void print_histogram(Histogram hist) {
    int32_t ibin, idim, totalbin;
    printf("mins:");
    for (idim = 0; idim < hist.ndims; ++idim)
        printf("%lf ", hist.mins[idim]);
    printf("   ");
    printf("maxs:");
    for (idim = 0; idim < hist.ndims; ++idim)
        printf("%lf ", hist.maxs[idim]);
    printf("   ");
    printf("percentinrange:%lf    ", hist.percentinrange);
    if (hist.issparse) {
        printf("sparse:%d|", hist.nnonemptybins);
        for (ibin = 0; ibin < hist.nnonemptybins; ++ibin) {
            printf("%d_%d,", hist.buffer[2 * ibin], hist.buffer[2 * ibin + 1]);
        }
    } else {
        printf("dense:");
        totalbin = total_number_of_bins(hist.ndims, hist.nbins);
        for (ibin = 0; ibin < totalbin; ++ibin) {
            printf("%d,", hist.buffer[ibin]);
        }
    }
    printf("\n");
}

void print_domain_histograms(
        int32_t ngridx, int32_t ngridy, int32_t ngridz,
        int32_t nhistx, int32_t nhisty, int32_t nhistz,
        int32_t ndims, int32_t nbins[], double log_base[],
        Histogram hists[], int32_t nhists) {
#ifdef PRINT_HISTOGRAM
    int32_t ihist;
    print_histogram_meta(
            ngridx, ngridy, ngridz, nhistx, nhisty, nhistz, ndims, nbins, log_base);
    for (ihist = 0; ihist < nhists; ++ihist)
        print_histogram(hists[ihist]);
#endif // PRINT_HISTOGRAM
}

void write_histogram_meta(
        FILE* outfile,
        int32_t ngridx, int32_t ngridy, int32_t ngridz,
        int32_t nhistx, int32_t nhisty, int32_t nhistz,
        int32_t ndims, int32_t nbins[], double log_base[]) {
    fwrite(&ndims, sizeof(int32_t), 1, outfile);
    fwrite(&ngridx, sizeof(int32_t), 1, outfile);
    fwrite(&ngridy, sizeof(int32_t), 1, outfile);
    fwrite(&ngridz, sizeof(int32_t), 1, outfile);
    fwrite(&nhistx, sizeof(int32_t), 1, outfile);
    fwrite(&nhisty, sizeof(int32_t), 1, outfile);
    fwrite(&nhistz, sizeof(int32_t), 1, outfile);
    fwrite(nbins, sizeof(int32_t), ndims, outfile);
    fwrite(log_base, sizeof(double), ndims, outfile);
}

void write_histogram(FILE* outfile, Histogram hist) {
    int32_t totalbin;
    fwrite(&hist.issparse, sizeof(int32_t), 1, outfile);
    fwrite(hist.mins, sizeof(double), hist.ndims, outfile);
    fwrite(hist.maxs, sizeof(double), hist.ndims, outfile);
    fwrite(&hist.percentinrange, sizeof(double), 1, outfile);
    fwrite(&hist.nnonemptybins, sizeof(int32_t), 1, outfile);
    if (hist.issparse) {
        fwrite(hist.buffer, sizeof(int32_t), 2 * hist.nnonemptybins, outfile);
    } else {
        totalbin = total_number_of_bins(hist.ndims, hist.nbins);
        fwrite(hist.buffer, sizeof(int32_t), totalbin, outfile);
    }
}

void write_domain_histograms(
        int32_t ngridx, int32_t ngridy, int32_t ngridz,
        int32_t nhistx, int32_t nhisty, int32_t nhistz,
        int32_t ndims, int32_t nbins[], double log_base[],
        Histogram hists[], int32_t nhists,
        double timestep, int32_t myid, int32_t iconfig) {
    int32_t ihist;
    char filename[100];
    FILE* outfile;
    // printf("../data/tracer-%.4E/pdfs-%03d.%05d\n", timestep, iconfig, myid);
    sprintf(filename, "../data/tracer-%.4E/pdfs-%03d.%05d", timestep, iconfig, myid);
    outfile = fopen(filename, "wb");
    write_histogram_meta(
            outfile, ngridx, ngridy, ngridz, nhistx, nhisty, nhistz, ndims, nbins, log_base);
    for (ihist = 0; ihist < nhists; ++ihist)
        write_histogram(outfile, hists[ihist]);
    fclose(outfile);
}

SamplingRegion construct_sampling_region(
        double* valuearrays[], int32_t ndims,
        int32_t ngridx, int32_t ngridy, int32_t ngridz,
        int32_t nhistx, int32_t nhisty, int32_t nhistz,
        int32_t ihistx, int32_t ihisty, int32_t ihistz) {
    int32_t i;
    int32_t ngridxperhist = ngridx / nhistx;
    int32_t ngridyperhist = ngridy / nhisty;
    int32_t ngridzperhist = ngridz / nhistz;
    SamplingRegion samplingregion;
    samplingregion.ndims = ndims;
    for (i = 0; i < ndims; ++i)
        samplingregion.valuearrays[i] = valuearrays[i];
    samplingregion.ngridx = ngridx;
    samplingregion.ngridy = ngridy;
    samplingregion.ngridz = ngridz;
    samplingregion.xbeg = (ihistx + 0) * ngridxperhist;
    samplingregion.xend = (ihistx + 1) * ngridxperhist;
    samplingregion.ybeg = (ihisty + 0) * ngridyperhist;
    samplingregion.yend = (ihisty + 1) * ngridyperhist;
    samplingregion.zbeg = (ihistz + 0) * ngridzperhist;
    samplingregion.zend = (ihistz + 1) * ngridzperhist;
    return samplingregion;
}

int32_t ids_to_flat(int32_t ndims, int32_t dims[], int32_t ids[]) {
    int32_t sum = 0, isum, prod, iprod;
    for (isum = 0; isum < ndims; ++isum) {
        prod = ids[isum];
        for (iprod = 0; iprod < isum; ++iprod)
            prod *= dims[iprod];
        sum += prod;
    }
    return sum;
}

int32_t values_to_bin_index(int32_t ndims, double values[],
        double mins[], double maxs[], int32_t nbins[]) {
    int32_t ibins[MAX_DIM], i;
    double scaled;
    for (i = 0; i < ndims; ++i) {
        if (values[i] < mins[i] || values[i] > maxs[i])
            return -1;
        scaled = (values[i] - mins[i]) / (maxs[i] - mins[i]);
        /// TODO: implement log scaling.
        ibins[i] = fmin(floor(scaled * nbins[i]), (double)(nbins[i] - 1));
    }
    // ibins to flat index
    return ids_to_flat(ndims, nbins, ibins);
}

int cmpdouble(const void *a, const void *b) {
    double aa = *(double*)a;
    double bb = *(double*)b;
    if (aa == bb) return 0;
    if (aa > bb) return 1;
    return -1;
}

void adjust_range_by_methods(int32_t ndims, SamplingRegion samplingregion,
        char* methods[], double mins[], double maxs[], double adjmins[],
        double adjmaxs[]) {
    // helper variables
    int32_t ngridx = samplingregion.ngridx;
    int32_t ngridy = samplingregion.ngridy;
    int32_t ngridz = samplingregion.ngridz;
    int32_t x, y, z;
    int32_t ibin, igrid, nnonemptybins, idim, isample, minpos, maxpos;
    double value, samplemin, samplemax, *samples;
    int32_t regionngrids =
            (samplingregion.xend - samplingregion.xbeg) *
            (samplingregion.yend - samplingregion.ybeg) *
            (samplingregion.zend - samplingregion.zbeg);
    for (idim = 0; idim < ndims; ++idim) {
        if (0 == strncmp(methods[idim], "range", 5)) {
            // do nothing
            adjmins[idim] = mins[idim];
            adjmaxs[idim] = maxs[idim];
        } else if (0 == strncmp(methods[idim], "normalized_range", 16)) {
            // get sample min and sample max by iterating all the samples
            samplemin = DBL_MAX;
            samplemax = -DBL_MAX;
            for (z = samplingregion.zbeg; z < samplingregion.zend; ++z)
            for (y = samplingregion.ybeg; y < samplingregion.yend; ++y)
            for (x = samplingregion.xbeg; x < samplingregion.xend; ++x) {
                igrid = x + y * ngridx + z * ngridy * ngridx;
                value = samplingregion.valuearrays[idim][igrid];
                samplemin = MIN(samplemin, value);
                samplemax = MAX(samplemax, value);
            }
            // linear interpolate to adjust the min and max
            adjmins[idim] = mins[idim] * (samplemax - samplemin) + samplemin;
            adjmaxs[idim] = maxs[idim] * (samplemax - samplemin) + samplemin;
        } else if (0 == strncmp(methods[idim], "percent_range", 13)) {
            // put all samples into an array for qsort
            samples = (double*)malloc(regionngrids * sizeof(double));
            isample = 0;
            for (z = samplingregion.zbeg; z < samplingregion.zend; ++z)
            for (y = samplingregion.ybeg; y < samplingregion.yend; ++y)
            for (x = samplingregion.xbeg; x < samplingregion.xend; ++x) {
                igrid = x + y * ngridx + z * ngridy * ngridx;
                samples[isample++] = samplingregion.valuearrays[idim][igrid];
            }
            // quicksort the array
            qsort(samples, regionngrids, sizeof(double), cmpdouble);
            // query the position in the sorted array for the min and max.
            minpos = floor(mins[idim] * regionngrids) + 0.5;
            maxpos = ceil(maxs[idim] * regionngrids) + 0.5;
            adjmins[idim] = samples[minpos];
            adjmaxs[idim] = samples[maxpos];
            // free the allocated array
            free(samples);
        }
    }
}

Histogram generate_histogram_from_sampling_region(
        int32_t ndims, SamplingRegion samplingregion, int32_t nbins[],
        double mins[], double maxs[]) {
    // helper variables
    int32_t ngridx = samplingregion.ngridx;
    int32_t ngridy = samplingregion.ngridy;
    int32_t ngridz = samplingregion.ngridz;
    int32_t x, y, z;
    int32_t ibin, igrid, nnonemptybins, idim;
    double values[MAX_DIM], scaled;
    int32_t regionngrids =
            (samplingregion.xend - samplingregion.xbeg) *
            (samplingregion.yend - samplingregion.ybeg) *
            (samplingregion.zend - samplingregion.zbeg);
    int32_t valuesinrange = regionngrids;
    double percentinrange;
    int32_t totalbin = total_number_of_bins(ndims, nbins);
    int32_t* frequencies = (int32_t*)calloc(totalbin, sizeof(int32_t));
    // sample each grid point
    for (z = samplingregion.zbeg; z < samplingregion.zend; ++z)
    for (y = samplingregion.ybeg; y < samplingregion.yend; ++y)
    for (x = samplingregion.xbeg; x < samplingregion.xend; ++x) {
        igrid = x + y * ngridx + z * ngridy * ngridx;
        for (idim = 0; idim < ndims; ++idim)
            values[idim] = samplingregion.valuearrays[idim][igrid];
        ibin = values_to_bin_index(ndims, values, mins, maxs, nbins);
        if (ibin < 0)
            --valuesinrange;
        else
            ++frequencies[ibin];
    }
    percentinrange = (double)valuesinrange / (double)regionngrids;
    // construct a histogram
    Histogram hist = frequencies_to_histogram(
            ndims, frequencies, nbins, mins, maxs, percentinrange);
    // clean up
    free(frequencies);
    return hist;
}

/// TODO: consider exposing this to the fortran interface
void generate_and_output_histogram(
        double* values[], double log_bases[], char* methods[],
        double mins[], double maxs[],
        double timestep,
        int32_t ndims,
        int32_t ngridx, int32_t ngridy, int32_t ngridz,
        int32_t nhistx, int32_t nhisty, int32_t nhistz,
        int32_t nbins[], int32_t myid, int32_t iconfig) {
    // declare variables
    int32_t nhists = nhistx * nhisty * nhistz;
    int32_t ihist, ihistx, ihisty, ihistz;
    SamplingRegion samplingregion;
    Histogram* hists = malloc(nhists * sizeof(Histogram));
    double adjmins[MAX_DIM], adjmaxs[MAX_DIM];
    // generate the histograms
    for (ihistz = 0; ihistz < nhistz; ++ihistz)
    for (ihisty = 0; ihisty < nhisty; ++ihisty)
    for (ihistx = 0; ihistx < nhistx; ++ihistx) {
        ihist = ihistx + ihisty * nhistx + ihistz * nhisty * nhistx;
        samplingregion = construct_sampling_region(
                values, ndims, ngridx, ngridy, ngridz,
                nhistx, nhisty, nhistz, ihistx, ihisty, ihistz);
        adjust_range_by_methods(
                ndims, samplingregion, methods, mins, maxs, adjmins, adjmaxs);
        hists[ihist] = generate_histogram_from_sampling_region(
                ndims, samplingregion, nbins, adjmins, adjmaxs);
    }
    // print to screen and write to file
    print_domain_histograms(
            ngridx, ngridy, ngridz, nhistx, nhisty, nhistz,
            ndims, nbins, log_bases, hists, nhists);
    write_domain_histograms(
            ngridx, ngridy, ngridz, nhistx, nhisty, nhistz,
            ndims, nbins, log_bases, hists, nhists, timestep,
            myid, iconfig);
    // deallocate histograms
    for (ihist = 0; ihist < nhists; ++ihist)
        deallocate_histogram(hists[ihist]);
    free(hists);
}

void c_generate_and_output_histogram_1d_(
        double *ptr_values,
        double *ptr_log_base,
        char *ptr_method,
        double *ptr_min, double *ptr_max,
        double *ptr_timestep,
        int32_t *ptr_ngridx, int32_t *ptr_ngridy, int32_t *ptr_ngridz,
        int32_t *ptr_nhistx, int32_t *ptr_nhisty, int32_t *ptr_nhistz,
        int32_t *ptr_nbin,
        int32_t *ptr_myid, int32_t *ptr_rootid, MPI_Fint *ptr_comm,
        int32_t *ptr_iconfig) {
    // convert from pointers to actual input arguments
    double *values = ptr_values, log_base = *ptr_log_base;
    char *method = ptr_method;
    double minimum = *ptr_min, maximum = *ptr_max;
    double timestep = *ptr_timestep;
    int32_t ngridx = *ptr_ngridx, ngridy = *ptr_ngridy, ngridz = *ptr_ngridz;
    int32_t nhistx = *ptr_nhistx, nhisty = *ptr_nhisty, nhistz = *ptr_nhistz;
    int32_t nbins = *ptr_nbin;
    int32_t myid = *ptr_myid, rootid = *ptr_rootid;
    MPI_Comm comm = MPI_Comm_f2c(*ptr_comm);
    int32_t iconfig = *ptr_iconfig;
    // generate and output
    generate_and_output_histogram(&values, &log_base, &method, &minimum,
            &maximum, timestep, 1, ngridx, ngridy, ngridz, nhistx, nhisty,
            nhistz, &nbins, myid, iconfig);
}

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
        int32_t *ptr_iconfig) {
    // convert from pointers to actual input arguments
    double* values[MAX_DIM] = {ptr_values_x, ptr_values_y};
    double log_bases[MAX_DIM] = {*ptr_log_base_x, *ptr_log_base_y};
    char* methods[MAX_DIM] = {ptr_method_x, ptr_method_y};
    double mins[MAX_DIM] = {*ptr_min_x, *ptr_min_y};
    double maxs[MAX_DIM] = {*ptr_max_x, *ptr_max_y};
    double timestep = *ptr_timestep;
    int32_t ngridx = *ptr_ngridx, ngridy = *ptr_ngridy, ngridz = *ptr_ngridz;
    int32_t nhistx = *ptr_nhistx, nhisty = *ptr_nhisty, nhistz = *ptr_nhistz;
    int32_t nbins[MAX_DIM] = {*ptr_nbinx, *ptr_nbiny};
    int32_t myid = *ptr_myid, rootid = *ptr_rootid;
    MPI_Comm comm = MPI_Comm_f2c(*ptr_comm);
    int32_t iconfig = *ptr_iconfig;
    // generate and output
    generate_and_output_histogram(values, log_bases, methods, mins, maxs,
            timestep, 2, ngridx, ngridy, ngridz, nhistx, nhisty, nhistz, nbins,
            myid, iconfig);
}

void c_generate_and_output_histogram_3d_(
        double *ptr_values_x, double *ptr_values_y, double *ptr_values_z,
        double *ptr_log_base_x, double *ptr_log_base_y, double *ptr_log_base_z,
        char *ptr_method_x, char *ptr_method_y, char *ptr_method_z,
        double *ptr_min_x, double *ptr_min_y, double *ptr_min_z,
        double *ptr_max_x, double *ptr_max_y, double *ptr_max_z,
        double *ptr_timestep,
        int32_t *ptr_ngridx, int32_t *ptr_ngridy, int32_t *ptr_ngridz,
        int32_t *ptr_nhistx, int32_t *ptr_nhisty, int32_t *ptr_nhistz,
        int32_t *ptr_nbinx, int32_t *ptr_nbiny, int32_t *ptr_nbinz,
        int32_t *ptr_myid, int32_t *ptr_rootid, MPI_Fint *ptr_comm,
        int32_t *ptr_iconfig) {
    // convert from pointers to actual input arguments
    double* values[MAX_DIM] = {ptr_values_x, ptr_values_y, ptr_values_z};
    double log_bases[MAX_DIM] =
            {*ptr_log_base_x, *ptr_log_base_y, *ptr_log_base_z};
    char* methods[MAX_DIM] = {ptr_method_x, ptr_method_y, ptr_method_z};
    double mins[MAX_DIM] = {*ptr_min_x, *ptr_min_y, *ptr_min_z};
    double maxs[MAX_DIM] = {*ptr_max_x, *ptr_max_y, *ptr_max_z};
    double timestep = *ptr_timestep;
    int32_t ngridx = *ptr_ngridx, ngridy = *ptr_ngridy, ngridz = *ptr_ngridz;
    int32_t nhistx = *ptr_nhistx, nhisty = *ptr_nhisty, nhistz = *ptr_nhistz;
    int32_t nbins[MAX_DIM] = {*ptr_nbinx, *ptr_nbiny, *ptr_nbinz};
    int32_t myid = *ptr_myid, rootid = *ptr_rootid;
    MPI_Comm comm = MPI_Comm_f2c(*ptr_comm);
    int32_t iconfig = *ptr_iconfig;
    // generate and output
    generate_and_output_histogram(values, log_bases, methods, mins, maxs,
            timestep, 3, ngridx, ngridy, ngridz, nhistx, nhisty, nhistz, nbins,
            myid, iconfig);
}
