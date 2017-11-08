#include "histogram.h"
#include "histogram_in.h"
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

static char dirpre[100] = "../data_pdf/pdf-";

/// TODO: change histogram dimension to a, b, c

// #define PRINT_HISTOGRAM

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
        int32_t ndims, int32_t* frequencies, int32_t nbins[],
        int32_t nnonemptybins, double mins[], double maxs[],
        double percentinrange) {
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
        int32_t ndims, int32_t* frequencies, int32_t nbins[],
        int32_t nnonemptybins, double mins[], double maxs[],
        double percentinrange) {
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
        return frequencies_to_histogram_dense(ndims, frequencies, nbins,
                nnonemptybins, mins, maxs, percentinrange);
    }
    return frequencies_to_histogram_sparse(ndims, frequencies, nbins,
            nnonemptybins, mins, maxs, percentinrange);
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
    print_histogram_meta(ngridx, ngridy, ngridz, nhistx, nhisty, nhistz, ndims,
            nbins, log_base);
    for (ihist = 0; ihist < nhists; ++ihist)
        print_histogram(hists[ihist]);
#endif // PRINT_HISTOGRAM
}

void write_histogram_meta(
        FILE* outfile,
        int32_t ngridx, int32_t ngridy, int32_t ngridz,
        int32_t nhistx, int32_t nhisty, int32_t nhistz,
        int32_t ndims, int32_t nbins[], double log_base[]) {
    char* buffer;
    int32_t nbyte;
    serialize_histogram_meta(ngridx, ngridy, ngridz, nhistx, nhisty, nhistz,
            ndims, log_base, &buffer, &nbyte);
    fwrite(buffer, sizeof(char), nbyte, outfile);
    free(buffer);

    // fwrite(&ndims, sizeof(int32_t), 1, outfile);
    // fwrite(&ngridx, sizeof(int32_t), 1, outfile);
    // fwrite(&ngridy, sizeof(int32_t), 1, outfile);
    // fwrite(&ngridz, sizeof(int32_t), 1, outfile);
    // fwrite(&nhistx, sizeof(int32_t), 1, outfile);
    // fwrite(&nhisty, sizeof(int32_t), 1, outfile);
    // fwrite(&nhistz, sizeof(int32_t), 1, outfile);
    // fwrite(log_base, sizeof(double), ndims, outfile);
}

void serialize_histogram_meta(
        int32_t ngridx, int32_t ngridy, int32_t ngridz,
        int32_t nhistx, int32_t nhisty, int32_t nhistz,
        int32_t ndims, double log_base[], char** buffer, int32_t* nbyte) {
    int32_t ibyte;
    *buffer = malloc(7 * sizeof(int32_t) + ndims * sizeof(double));
    *nbyte = 7 * sizeof(int32_t) + ndims * sizeof(double);
    ibyte = 0;
    memcpy(*buffer, &ndims, sizeof(int32_t)); ibyte += sizeof(int32_t);
    memcpy(*buffer + ibyte, &ngridx, sizeof(int32_t)); ibyte += sizeof(int32_t);
    memcpy(*buffer + ibyte, &ngridy, sizeof(int32_t)); ibyte += sizeof(int32_t);
    memcpy(*buffer + ibyte, &ngridz, sizeof(int32_t)); ibyte += sizeof(int32_t);
    memcpy(*buffer + ibyte, &nhistx, sizeof(int32_t)); ibyte += sizeof(int32_t);
    memcpy(*buffer + ibyte, &nhisty, sizeof(int32_t)); ibyte += sizeof(int32_t);
    memcpy(*buffer + ibyte, &nhistz, sizeof(int32_t)); ibyte += sizeof(int32_t);
    memcpy(*buffer + ibyte, log_base, sizeof(double) * ndims);
}

void deserialize_histogram_meta(
        char* buffer, int32_t nbyte,
        int32_t* ngridx, int32_t* ngridy, int32_t* ngridz,
        int32_t* nhistx, int32_t* nhisty, int32_t* nhistz,
        int32_t* ndims, double log_base[]) {
    int32_t ibyte;
    ibyte = 0;
    memcpy(ndims, buffer + ibyte, sizeof(int32_t)); ibyte += sizeof(int32_t);
    memcpy(ngridx, buffer + ibyte, sizeof(int32_t)); ibyte += sizeof(int32_t);
    memcpy(ngridy, buffer + ibyte, sizeof(int32_t)); ibyte += sizeof(int32_t);
    memcpy(ngridz, buffer + ibyte, sizeof(int32_t)); ibyte += sizeof(int32_t);
    memcpy(nhistx, buffer + ibyte, sizeof(int32_t)); ibyte += sizeof(int32_t);
    memcpy(nhisty, buffer + ibyte, sizeof(int32_t)); ibyte += sizeof(int32_t);
    memcpy(nhistz, buffer + ibyte, sizeof(int32_t)); ibyte += sizeof(int32_t);
    memcpy(log_base, buffer + ibyte, sizeof(double) * (*ndims));
}

void write_histogram(FILE* outfile, Histogram hist) {
    char* buffer;
    int32_t nbyte;
    serialize_histogram(hist, &buffer, &nbyte);
    fwrite(buffer, sizeof(char), nbyte, outfile);
    free(buffer);

    // int32_t totalbin;
    // fwrite(&hist.issparse, sizeof(int32_t), 1, outfile);
    // fwrite(hist.mins, sizeof(double), hist.ndims, outfile);
    // fwrite(hist.maxs, sizeof(double), hist.ndims, outfile);
    // fwrite(hist.nbins, sizeof(int32_t), hist.ndims, outfile);
    // fwrite(&hist.percentinrange, sizeof(double), 1, outfile);
    // fwrite(&hist.nnonemptybins, sizeof(int32_t), 1, outfile);
    // if (hist.issparse) {
    //     fwrite(hist.buffer, sizeof(int32_t), 2 * hist.nnonemptybins,
    //             outfile);
    // } else {
    //     totalbin = total_number_of_bins(hist.ndims, hist.nbins);
    //     fwrite(hist.buffer, sizeof(int32_t), totalbin, outfile);
    // }
}

int hist_buffer_byte_count(int32_t issparse, int32_t ndims, int32_t nbins[],
        int32_t nnonemptybins) {
    if (issparse) {
        return sizeof(int32_t) * 2 * nnonemptybins;
    }
    return sizeof(int32_t) * total_number_of_bins(ndims, nbins);
}

int hist_byte_count(int32_t issparse, int32_t ndims, int32_t nbins[],
        int32_t nnonemptybins) {
    return 2 * sizeof(int32_t) +
            2 * sizeof(double) * ndims +
            sizeof(int32_t) * ndims +
            sizeof(double) +
            hist_buffer_byte_count(issparse, ndims, nbins, nnonemptybins);
}

int buffer_to_hist_byte(char* buffer, int ndims) {
    int32_t issparse, *nbins, nnonemptybins, totalbin, idim, nBytes;
    issparse = *buffer;
    nbins = malloc(ndims * sizeof(int32_t));
    for (idim = 0; idim < ndims; ++idim) {
        nbins[idim] =
                *(buffer + (1 + idim) * sizeof(int32_t)
                    + 2 * ndims * sizeof(double));
    }
    nnonemptybins =
            *(buffer + (1 + ndims) * sizeof(int32_t)
                + (2 * ndims + 1) * sizeof(double));
    nBytes = hist_byte_count(issparse, ndims, nbins, nnonemptybins);
    free(nbins);
    return nBytes;
}

int buffer_to_hist_count(int nBytes, char* buffer, int ndims) {
    int offset, nHists;
    nHists = 0;
    offset = 0;
    while (offset < nBytes) {
        offset += buffer_to_hist_byte(buffer + offset, ndims);
        ++nHists;
    }
    return nHists;
}

void serialize_histogram(Histogram hist, char** buffer, int32_t *nbyte) {
    int32_t totalbin, ibyte;
    totalbin = total_number_of_bins(hist.ndims, hist.nbins);
    *nbyte =
            hist_byte_count(
                hist.issparse, hist.ndims, hist.nbins, hist.nnonemptybins);
    *buffer = malloc(*nbyte);
    ibyte = 0;
    memcpy(*buffer, &hist.issparse, sizeof(int32_t)); ibyte += sizeof(int32_t);
    memcpy(*buffer + ibyte, hist.mins, sizeof(double) * hist.ndims);
    ibyte += sizeof(double) * hist.ndims;
    memcpy(*buffer + ibyte, hist.maxs, sizeof(double) * hist.ndims);
    ibyte += sizeof(double) * hist.ndims;
    memcpy(*buffer + ibyte, hist.nbins, sizeof(int32_t) * hist.ndims);
    ibyte += sizeof(int32_t) * hist.ndims;
    memcpy(*buffer + ibyte, &hist.percentinrange, sizeof(double));
    ibyte += sizeof(double);
    memcpy(*buffer + ibyte, &hist.nnonemptybins, sizeof(int32_t));
    ibyte += sizeof(int32_t);
    memcpy(*buffer + ibyte, hist.buffer,
            hist_buffer_byte_count(
                hist.issparse, hist.ndims, hist.nbins, hist.nnonemptybins));
}

void serialize_histograms(
        int nHists, Histogram hists[], int* nBytes, char** buffer) {
    int32_t iHist, *nBytesPerHists, offset;
    char **buffers;
    nBytesPerHists = malloc(nHists * sizeof(int32_t));
    buffers = malloc(nHists * sizeof(char*));
    for (iHist = 0; iHist < nHists; ++iHist) {
        serialize_histogram(
                hists[iHist], &buffers[iHist], &nBytesPerHists[iHist]);
    }
    *nBytes = 0;
    for (iHist = 0; iHist < nHists; ++iHist) {
        *nBytes += nBytesPerHists[iHist];
    }
    *buffer = malloc((*nBytes) * sizeof(char));
    offset = 0;
    for (iHist = 0; iHist < nHists; ++iHist) {
        memcpy(*buffer + offset, buffers[iHist], nBytesPerHists[iHist]);
        offset += nBytesPerHists[iHist];
    }
    for (iHist = 0; iHist < nHists; ++iHist) {
        free(buffers[iHist]);
    }
    free(nBytesPerHists);
    free(buffers);
}

void deserialize_histogram(
        char* buffer, int ndims, int* histByteCount, Histogram* hist) {
    int offset, histBufferByteCount; 
    offset = 0;
    hist->ndims = ndims;
    memcpy(&hist->issparse, buffer, sizeof(int32_t)); offset += sizeof(int32_t);
    memcpy(hist->mins, buffer + offset, ndims * sizeof(double));
    offset += ndims * sizeof(double);
    memcpy(hist->maxs, buffer + offset, ndims * sizeof(double));
    offset += ndims * sizeof(double);
    memcpy(hist->nbins, buffer + offset, ndims * sizeof(int32_t));
    offset += ndims * sizeof(int32_t);
    memcpy(&hist->percentinrange, buffer + offset, sizeof(double));
    offset += sizeof(double);
    memcpy(&hist->nnonemptybins, buffer + offset, sizeof(int32_t));
    offset += sizeof(int32_t);
    histBufferByteCount =
            hist_buffer_byte_count(
                hist->issparse, hist->ndims, hist->nbins, hist->nnonemptybins);
    *histByteCount =
            hist_byte_count(
                hist->issparse, hist->ndims, hist->nbins, hist->nnonemptybins);
    hist->buffer = malloc(histBufferByteCount);
    memcpy(hist->buffer, buffer + offset, histBufferByteCount);
}

void deserialize_histograms(
        int nBytes, char* buffer, int ndims, int* nHists, Histogram** hists) {
    int iHist, offset, histBufferBytes;
    *nHists = buffer_to_hist_count(nBytes, buffer, ndims);
    *hists = malloc(*nHists * sizeof(Histogram));
    offset = 0;
    for (iHist = 0; iHist < *nHists; ++iHist) {
        deserialize_histogram(
                buffer + offset, ndims, &histBufferBytes, &(*hists)[iHist]);
        offset += histBufferBytes;
    }
}

int is_same_histogram(Histogram a, Histogram b) {
    int idim, ibyte, bufferByteCount;
    if (a.ndims != b.ndims) return 0;
    if (a.issparse != b.issparse) return 0;
    if (a.nnonemptybins != b.nnonemptybins) return 0;
    if (fabs(a.percentinrange - b.percentinrange) > 0.0001) return 0;
    for (idim = 0; idim < a.ndims; ++idim) {
        if (a.nbins[idim] != b.nbins[idim]) return 0;
        if (fabs(a.mins[idim] - b.mins[idim]) > 0.0001) return 0;
        if (fabs(a.maxs[idim] - b.maxs[idim]) > 0.0001) return 0;
    }
    bufferByteCount =
            hist_buffer_byte_count(
                a.issparse, a.ndims, a.nbins, a.nnonemptybins);
    for (ibyte = 0; ibyte < bufferByteCount; ++ibyte) {
        if (((char*)a.buffer)[ibyte] != ((char*)b.buffer)[ibyte]) return 0;
    }
    return 1;
}

void write_domain_histograms(
        int32_t ngridx, int32_t ngridy, int32_t ngridz,
        int32_t nhistx, int32_t nhisty, int32_t nhistz,
        int32_t ndims, int32_t nbins[], double log_base[],
        Histogram hists[], int32_t nhists,
        double timestep, int32_t yid, int32_t rootid, MPI_Comm comm,
        int32_t xzid, int32_t iconfig) {
    char* buffer;
    int32_t nbyte, ihist;
    char filename[100];
    FILE* outfile;
    sprintf(filename, "%s%.4E/pdfs-%03d.%05d", dirpre, timestep, iconfig, yid);
    // sprintf(filename, "../data/tracer-%.4E/pdfs-%03d.%05d", timestep,
            // iconfig, yid);
    outfile = fopen(filename, "wb");
    serialize_domain_histograms(ngridx, ngridy, ngridz, nhistx, nhisty, nhistz,
            ndims, log_base, hists, nhists, &buffer, &nbyte);
    fwrite(buffer, sizeof(char), nbyte, outfile);
    free(buffer);
    fclose(outfile);
}

void serialize_domain_histograms(
        int32_t ngridx, int32_t ngridy, int32_t ngridz,
        int32_t nhistx, int32_t nhisty, int32_t nhistz,
        int32_t ndims, double log_base[],
        Histogram hists[], int32_t nhists, char** buffer, int32_t* nbyte) {
    char *metabuffer, *histbuffer;
    int32_t metanbyte, histnbyte, ihist;
    serialize_histogram_meta(ngridx, ngridy, ngridz, nhistx, nhisty, nhistz,
            ndims, log_base, &metabuffer, &metanbyte);
    serialize_histograms(nhists, hists, &histnbyte, &histbuffer);
    *nbyte = (metanbyte + histnbyte) * sizeof(char);
    *buffer = malloc(*nbyte);
    memcpy(*buffer, metabuffer, metanbyte);
    memcpy(*buffer + metanbyte, histbuffer, histnbyte);
    // free memory
    free(metabuffer);
    free(histbuffer);
}

void gather_buffers(int nBytesPerRank, char* bufferPerRank,
        int** nBytesPerRanks, int* nBytes, char** buffer, int root,
        MPI_Comm comm) {
    int rank, nRanks, *displs, iRank;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nRanks);
    if (root == rank) {
        *nBytesPerRanks = malloc(nRanks * sizeof(int));
    }
    MPI_Gather(&nBytesPerRank, 1, MPI_INT, *nBytesPerRanks, 1, MPI_INT, root,
            comm);
    if (root == rank) {
        displs = malloc(nRanks * sizeof(int));
        *nBytes = 0;
        for (iRank = 0; iRank < nRanks; ++iRank) {
            displs[iRank] = *nBytes;
            *nBytes += (*nBytesPerRanks)[iRank];
        }
        *buffer = malloc(*nBytes);
    }
    MPI_Gatherv(bufferPerRank, nBytesPerRank, MPI_CHAR, *buffer,
            *nBytesPerRanks, displs, MPI_CHAR, root, comm);
    if (root == rank) {
        free(displs);
    }
}

void write_histograms_grouped_by_mpi_comm(
        int32_t ngridx, int32_t ngridy, int32_t ngridz,
        int32_t nhistx, int32_t nhisty, int32_t nhistz,
        int32_t ndims, double log_base[], Histogram hists[], int32_t nhists,
        double timestep, int32_t yid, int32_t rootid, int32_t ypes,
        MPI_Comm comm, int32_t xzid, int32_t iconfig, const char* midname) {
    int32_t ipes, domainnbyte, *domainnbytes, globalnbyte = 0;
    char filename[100], *domainbuffer, *globalbuffer;
    FILE* outfile;
    // serialize domain histograms
    serialize_domain_histograms(ngridx, ngridy, ngridz, nhistx, nhisty, nhistz,
            ndims, log_base, hists, nhists, &domainbuffer, &domainnbyte);
    // merge histograms in the y column
    gather_buffers(domainnbyte, domainbuffer, &domainnbytes, &globalnbyte,
            &globalbuffer, rootid, comm);
    // write to file
    if (yid == rootid) {
        sprintf(filename, "%s%.4E/pdfs-%s-%03d.%05d", dirpre, timestep, midname,
                iconfig, xzid);
        outfile = fopen(filename, "wb");
        fwrite((void*)globalbuffer, sizeof(char), globalnbyte, outfile);
        fclose(outfile);
    }
    // clean up memory
    if (yid == rootid) {
        free(domainnbytes);
        free(globalbuffer);
    }
    free(domainbuffer);
}

void write_ycolumn_histograms(
        int32_t ngridx, int32_t ngridy, int32_t ngridz,
        int32_t nhistx, int32_t nhisty, int32_t nhistz,
        int32_t ndims, double log_base[], Histogram hists[], int32_t nhists,
        double timestep, int32_t yid, int32_t rootid, int32_t ypes,
        MPI_Comm comm, int32_t xzid, int32_t iconfig) {
    write_histograms_grouped_by_mpi_comm(
            ngridx, ngridy, ngridz, nhistx, nhisty, nhistz, ndims, log_base,
            hists, nhists, timestep, yid, rootid, ypes, comm, xzid, iconfig,
            "ycolumn");
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

void flatten_sampling_region(
        SamplingRegion samplingregion, int32_t idim, double samples[]) {
    int32_t isample = 0, x, y, z, igrid;
    int32_t ngridx = samplingregion.ngridx;
    int32_t ngridy = samplingregion.ngridy;
    int32_t ngridz = samplingregion.ngridz;
    for (z = samplingregion.zbeg; z < samplingregion.zend; ++z)
    for (y = samplingregion.ybeg; y < samplingregion.yend; ++y)
    for (x = samplingregion.xbeg; x < samplingregion.xend; ++x) {
        igrid = x + y * ngridx + z * ngridy * ngridx;
        samples[isample++] = samplingregion.valuearrays[idim][igrid];
    }
}

int32_t flatten_sampling_region_with_range(SamplingRegion samplingregion,
        int32_t idim, double min, double max, double samples[]) {
    int32_t isample = 0, x, y, z, igrid;
    int32_t ngridx = samplingregion.ngridx;
    int32_t ngridy = samplingregion.ngridy;
    int32_t ngridz = samplingregion.ngridz;
    double value;
    for (z = samplingregion.zbeg; z < samplingregion.zend; ++z)
    for (y = samplingregion.ybeg; y < samplingregion.yend; ++y)
    for (x = samplingregion.xbeg; x < samplingregion.xend; ++x) {
        igrid = x + y * ngridx + z * ngridy * ngridx;
        value = samplingregion.valuearrays[idim][igrid];
        if (min <= value && value < max)
            samples[isample++] = value;
    }
    return isample;
}

int32_t count_samples_in_sampling_region_in_range(SamplingRegion samplingregion,
        int32_t idim, double min, double max) {
    int32_t isample = 0, x, y, z, igrid;
    int32_t ngridx = samplingregion.ngridx;
    int32_t ngridy = samplingregion.ngridy;
    int32_t ngridz = samplingregion.ngridz;
    double value;
    for (z = samplingregion.zbeg; z < samplingregion.zend; ++z)
    for (y = samplingregion.ybeg; y < samplingregion.yend; ++y)
    for (x = samplingregion.xbeg; x < samplingregion.xend; ++x) {
        igrid = x + y * ngridx + z * ngridy * ngridx;
        value = samplingregion.valuearrays[idim][igrid];
        if (min <= value && value < max)
            ++isample;
    }
    return isample;
}

void adjust_range_by_methods(int32_t ndims, SamplingRegion samplingregion,
        char* methods[], double mins[], double maxs[], double adjmins[],
        double adjmaxs[]) {
    // helper variables
    int32_t ngridx = samplingregion.ngridx;
    int32_t ngridy = samplingregion.ngridy;
    int32_t ngridz = samplingregion.ngridz;
    int32_t x, y, z, i;
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
            flatten_sampling_region(samplingregion, idim, samples);
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

void compute_nbins(int32_t ndims, SamplingRegion samplingregion, double mins[],
        double maxs[], char* nbinstrs[], int32_t nbins[]) {
    int32_t idim, regionngrids, i, nsamples;
    double *samples, iqr, q1, q3, h;
    for (idim = 0; idim < ndims; ++idim) {
        int32_t nbin = atoi(nbinstrs[idim]);
        if (nbin > 0) {
            nbins[idim] = nbin;
        } else if (0 == strncmp(nbinstrs[idim], "sturges", 7)) {
            nsamples = count_samples_in_sampling_region_in_range(
                    samplingregion, idim, mins[idim], maxs[idim]);
            // printf("nsamples = %d\n", nsamples);
            if (0 == nsamples) {
                nbins[idim] = 1;
            } else {
                nbins[idim] = ceil(log2((double)nsamples)) + 1;
            }
        } else { // if (0 == strncmp(nbinstrs[idim], "freedman", 8)) {
            // put all samples into an array for qsort
            regionngrids =
                    (samplingregion.xend - samplingregion.xbeg) *
                    (samplingregion.yend - samplingregion.ybeg) *
                    (samplingregion.zend - samplingregion.zbeg);
            samples = (double*)malloc(regionngrids * sizeof(double));
            nsamples = flatten_sampling_region_with_range(
                    samplingregion, idim, mins[idim], maxs[idim], samples);
            // quicksort the array
            qsort(samples, nsamples, sizeof(double), cmpdouble);
            // compute IQR (interquartile range)
            q1 = samples[(int)(1.0 / 4.0 * nsamples)];
            q3 = samples[(int)(3.0 / 4.0 * nsamples)];
            iqr = q3 - q1;
            if (0 == nsamples) {
                h = 0.0;
            } else {
                h = 2 * iqr / pow((double)nsamples, 1.0/3.0);
            }
            if (h < DBL_EPSILON) {
                nbins[idim] = 1;
            } else {
                nbins[idim] = MIN(MAX_NBINS, (maxs[idim] - mins[idim]) / h);
            }
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
        char* nbinstrs[], int32_t yid, int32_t rootid, int32_t ypes,
        MPI_Comm comm, int32_t xzid, int32_t iconfig) {
    // declare variables
    int32_t nhists = nhistx * nhisty * nhistz;
    int32_t ihist, ihistx, ihisty, ihistz;
    SamplingRegion samplingregion;
    Histogram* hists = malloc(nhists * sizeof(Histogram));
    double adjmins[MAX_DIM], adjmaxs[MAX_DIM];
    int32_t nbins[MAX_DIM];
    // generate the histograms
    for (ihistz = 0; ihistz < nhistz; ++ihistz)
    for (ihisty = 0; ihisty < nhisty; ++ihisty)
    for (ihistx = 0; ihistx < nhistx; ++ihistx) {
        ihist = ihistx + ihisty * nhistx + ihistz * nhisty * nhistx;
        samplingregion = construct_sampling_region(
                values, ndims, ngridx, ngridy, ngridz,
                nhistx, nhisty, nhistz, ihistx, ihisty, ihistz);
        /// TODO: optimize the following 3 functions into first flattening the
        /// sampling region into an array, and then create histogram with the 1d
        /// array.
        adjust_range_by_methods(
                ndims, samplingregion, methods, mins, maxs, adjmins, adjmaxs);
        compute_nbins(ndims, samplingregion, adjmins, adjmaxs, nbinstrs, nbins);
        hists[ihist] = generate_histogram_from_sampling_region(
                ndims, samplingregion, nbins, adjmins, adjmaxs);
    }
    // print to screen and write to file
    print_domain_histograms(
            ngridx, ngridy, ngridz, nhistx, nhisty, nhistz,
            ndims, nbins, log_bases, hists, nhists);
    write_ycolumn_histograms(
            ngridx, ngridy, ngridz, nhistx, nhisty, nhistz,
            ndims, log_bases, hists, nhists, timestep,
            yid, rootid, ypes, comm, xzid, iconfig);
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
        char *ptr_nbin,
        int32_t *ptr_yid, int32_t *ptr_rootid, int32_t *ptr_ypes,
        MPI_Fint *ptr_comm, int32_t *ptr_xzid, int32_t *ptr_iconfig) {
    // convert from pointers to actual input arguments
    double *values = ptr_values, log_base = *ptr_log_base;
    char *method = ptr_method;
    double minimum = *ptr_min, maximum = *ptr_max;
    double timestep = *ptr_timestep;
    int32_t ngridx = *ptr_ngridx, ngridy = *ptr_ngridy, ngridz = *ptr_ngridz;
    int32_t nhistx = *ptr_nhistx, nhisty = *ptr_nhisty, nhistz = *ptr_nhistz;
    char *nbins = ptr_nbin;
    int32_t yid = *ptr_yid, rootid = *ptr_rootid, ypes = *ptr_ypes;
    MPI_Comm comm = MPI_Comm_f2c(*ptr_comm);
    int32_t xzid = *ptr_xzid;
    int32_t iconfig = *ptr_iconfig;
    // generate and output
    generate_and_output_histogram(&values, &log_base, &method, &minimum,
            &maximum, timestep, 1, ngridx, ngridy, ngridz, nhistx, nhisty,
            nhistz, &nbins, yid, rootid, ypes, comm, xzid, iconfig);
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
        char *ptr_nbinx, char *ptr_nbiny,
        int32_t *ptr_yid, int32_t *ptr_rootid, int32_t *ptr_ypes,
        MPI_Fint *ptr_comm, int32_t *ptr_xzid, int32_t *ptr_iconfig) {
    // convert from pointers to actual input arguments
    double* values[MAX_DIM] = {ptr_values_x, ptr_values_y};
    double log_bases[MAX_DIM] = {*ptr_log_base_x, *ptr_log_base_y};
    char* methods[MAX_DIM] = {ptr_method_x, ptr_method_y};
    double mins[MAX_DIM] = {*ptr_min_x, *ptr_min_y};
    double maxs[MAX_DIM] = {*ptr_max_x, *ptr_max_y};
    double timestep = *ptr_timestep;
    int32_t ngridx = *ptr_ngridx, ngridy = *ptr_ngridy, ngridz = *ptr_ngridz;
    int32_t nhistx = *ptr_nhistx, nhisty = *ptr_nhisty, nhistz = *ptr_nhistz;
    char* nbins[MAX_DIM] = {ptr_nbinx, ptr_nbiny};
    int32_t yid = *ptr_yid, rootid = *ptr_rootid, ypes = *ptr_ypes;
    MPI_Comm comm = MPI_Comm_f2c(*ptr_comm);
    int32_t xzid = *ptr_xzid;
    int32_t iconfig = *ptr_iconfig;
    // generate and output
    generate_and_output_histogram(values, log_bases, methods, mins, maxs,
            timestep, 2, ngridx, ngridy, ngridz, nhistx, nhisty, nhistz, nbins,
            yid, rootid, ypes, comm, xzid, iconfig);
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
        char *ptr_nbinx, char *ptr_nbiny, char *ptr_nbinz,
        int32_t *ptr_yid, int32_t *ptr_rootid, int32_t *ptr_ypes,
        MPI_Fint *ptr_comm, int32_t *ptr_xzid, int32_t *ptr_iconfig) {
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
    char* nbins[MAX_DIM] = {ptr_nbinx, ptr_nbiny, ptr_nbinz};
    int32_t yid = *ptr_yid, rootid = *ptr_rootid, ypes = *ptr_ypes;
    MPI_Comm comm = MPI_Comm_f2c(*ptr_comm);
    int32_t xzid = *ptr_xzid;
    int32_t iconfig = *ptr_iconfig;
    // generate and output
    generate_and_output_histogram(values, log_bases, methods, mins, maxs,
            timestep, 3, ngridx, ngridy, ngridz, nhistx, nhisty, nhistz, nbins,
            yid, rootid, ypes, comm, xzid, iconfig);
}

void c_generate_and_output_histogram_(
        double *ptr_values, double *ptr_log_bases, char* ptr_methods,
        char *ptr_nbins, double *ptr_mins, double *ptr_maxs,
        double *ptr_timestep, int32_t *ptr_ndims,
        int32_t *ptr_ngridx, int32_t *ptr_ngridy, int32_t *ptr_ngridz,
        int32_t *ptr_nhistx, int32_t *ptr_nhisty, int32_t *ptr_nhistz,
        int32_t *ptr_yid, int32_t *ptr_rootid,
        int32_t *ptr_ypes, MPI_Fint *ptr_comm, int32_t *ptr_xzid,
        int32_t *ptr_iconfig) {
    // constant variables
    int32_t stringLength = 30;
    // temporary variables
    int i;
    // convert from pointers to actual input arguments
    int32_t ndims = *ptr_ndims;
    double* log_bases = ptr_log_bases;
    char* methods[MAX_DIM];
    char* nbins[MAX_DIM];
    for (i = 0; i < ndims; ++i) {
        methods[i] = ptr_methods + stringLength * i;
        nbins[i] = ptr_nbins + stringLength + i;
    }
    double* mins = ptr_mins;
    double* maxs = ptr_maxs;
    double timestep = *ptr_timestep;
    int32_t ngridx = *ptr_ngridx, ngridy = *ptr_ngridy, ngridz = *ptr_ngridz;
    int32_t nhistx = *ptr_nhistx, nhisty = *ptr_nhisty, nhistz = *ptr_nhistz;
    int32_t yid = *ptr_yid, rootid = *ptr_rootid, ypes = *ptr_ypes;
    MPI_Comm comm = MPI_Comm_f2c(*ptr_comm);
    int32_t xzid = *ptr_xzid;
    int32_t iconfig = *ptr_iconfig;
    // values
    int32_t ngrid = ngridx * ngridy * ngridz;
    double* values[MAX_DIM];
    for (i = 0; i < ndims; ++i) {
        values[i] = ptr_values + ngrid * i;
    }
    // generate and output
    generate_and_output_histogram(values, log_bases, methods, mins, maxs,
            timestep, ndims, ngridx, ngridy, ngridz, nhistx, nhisty, nhistz,
            nbins, yid, rootid, ypes, comm, xzid, iconfig);
}