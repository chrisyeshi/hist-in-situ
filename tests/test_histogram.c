#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <histogram.h>
#include <histogram_in.h>

void test_total_number_of_bins() {
    int32_t d1[] = {4};
    int32_t d2[] = {3, 4};
    int32_t d3[] = {3, 4, 5};
    assert(4 == total_number_of_bins(1, d1));
    assert(12 == total_number_of_bins(2, d2));
    assert(60 == total_number_of_bins(3, d3));
}

void test_get_number_of_nonempty_bins() {
    int32_t three[5] = {0, 0, 1, 2, 3};
    assert(3 == get_number_of_nonempty_bins(three, 5));
    int32_t zero[5] = {0, 0, 0, 0, 0};
    assert(0 == get_number_of_nonempty_bins(zero, 5));
    assert(0 == get_number_of_nonempty_bins(0, 0));
}

void test_frequencies_to_histogram_1d_dense() {
    int32_t i;
    int32_t frequencies[5] = {0, 0, 1, 2, 3};
    int32_t nbins[] = {5};
    double mins[] = {12.34};
    double maxs[] = {34.56};
    Histogram hist = frequencies_to_histogram_dense(
            1, frequencies, nbins, 3, mins, maxs, 0.78);
    assert(1 == hist.ndims);
    assert(0 == hist.issparse);
    assert(5 == hist.nbins[0]);
    assert(3 == hist.nnonemptybins);
    assert(12.34 == hist.mins[0]);
    assert(34.56 == hist.maxs[0]);
    assert(0.78 == hist.percentinrange);
    assert(frequencies != hist.buffer);
    for (i = 0; i < 5; ++i)
        assert(hist.buffer[i] == frequencies[i]);
}

void test_frequencies_to_histogram_2d_dense() {
    int32_t i;
    int32_t frequencies[6] = {3, 0, 4, 0, 2, 45};
    int32_t nbins[] = {2, 3};
    double mins[] = {12.34, 23.45};
    double maxs[] = {34.56, 45.67};
    Histogram hist = frequencies_to_histogram_dense(
            2, frequencies, nbins, 4, mins, maxs, 0.67);
    assert(2 == hist.ndims);
    assert(0 == hist.issparse);
    assert(2 == hist.nbins[0]);
    assert(3 == hist.nbins[1]);
    assert(12.34 == hist.mins[0]);
    assert(23.45 == hist.mins[1]);
    assert(34.56 == hist.maxs[0]);
    assert(45.67 == hist.maxs[1]);
    assert(0.67 == hist.percentinrange);
    assert(frequencies != hist.buffer);
    for (i = 0; i < 6; ++i)
        assert(hist.buffer[i] == frequencies[i]);
}

void test_frequencies_to_histogram_3d_dense() {
    int32_t i;
    int32_t frequencies[8] = {0, 1, 2, 3, 4, 5, 6, 7};
    int32_t nbins[] = {2, 2, 2};
    double mins[] = {2.3, 3.4, 4.5};
    double maxs[] = {3.4, 4.5, 5.6};
    Histogram hist = frequencies_to_histogram_dense(
            3, frequencies, nbins, 7, mins, maxs, 0.56);
    assert(3 == hist.ndims);
    assert(0 == hist.issparse);
    assert(2 == hist.nbins[0]);
    assert(2 == hist.nbins[1]);
    assert(2 == hist.nbins[2]);
    assert(2.3 == hist.mins[0]);
    assert(3.4 == hist.mins[1]);
    assert(4.5 == hist.mins[2]);
    assert(3.4 == hist.maxs[0]);
    assert(4.5 == hist.maxs[1]);
    assert(5.6 == hist.maxs[2]);
    assert(0.56 == hist.percentinrange);
    assert(frequencies != hist.buffer);
    for (i = 0; i < 8; ++i)
        assert(hist.buffer[i] == frequencies[i]);
}

void test_frequencies_to_histogram_1d_sparse() {
    int32_t frequencies[7] = {3, 0, 0, 45, 0, 8, 0};
    int32_t nbins[] = {7};
    double mins[] = {5.6};
    double maxs[] = {6.7};
    Histogram hist = frequencies_to_histogram_sparse(
            1, frequencies, nbins, 3, mins, maxs, 0.56);
    assert(1 == hist.ndims);
    assert(1 == hist.issparse);
    assert(7 == hist.nbins[0]);
    assert(3 == hist.nnonemptybins);
    assert(5.6 == hist.mins[0]);
    assert(6.7 == hist.maxs[0]);
    assert(0.56 == hist.percentinrange);
    assert(frequencies != hist.buffer);
    assert(0 == hist.buffer[0]);
    assert(3 == hist.buffer[1]);
    assert(3 == hist.buffer[2]);
    assert(45 == hist.buffer[3]);
    assert(5 == hist.buffer[4]);
    assert(8 == hist.buffer[5]);
}

void test_frequencies_to_histogram_2d_sparse() {
    int32_t frequencies[4] = {0, 5, 0, 0};
    int32_t nbins[] = {2, 2};
    double mins[] = {12.34, 23.45};
    double maxs[] = {34.56, 45.67};
    Histogram hist = frequencies_to_histogram_sparse(
            2, frequencies, nbins, 1, mins, maxs, 0.78);
    assert(2 == hist.ndims);
    assert(1 == hist.issparse);
    assert(2 == hist.nbins[0]);
    assert(2 == hist.nbins[1]);
    assert(1 == hist.nnonemptybins);
    assert(12.34 == hist.mins[0]);
    assert(23.45 == hist.mins[1]);
    assert(34.56 == hist.maxs[0]);
    assert(45.67 == hist.maxs[1]);
    assert(0.78 == hist.percentinrange);
    assert(frequencies != hist.buffer);
    assert(1 == hist.buffer[0]);
    assert(5 == hist.buffer[1]);
}

void test_frequencies_to_histogram_3d_sparse() {
    int32_t frequencies[8] = {0, 9, 0, 0, 2, 0, 50, 0};
    int32_t nbins[] = {2, 2, 2};
    double mins[] = {1.2, 2.3, 3.4};
    double maxs[] = {2.3, 3.4, 4.5};
    Histogram hist = frequencies_to_histogram_sparse(
            3, frequencies, nbins, 3, mins, maxs, 0.45);
    assert(3 == hist.ndims);
    assert(1 == hist.issparse);
    assert(2 == hist.nbins[0]);
    assert(2 == hist.nbins[1]);
    assert(2 == hist.nbins[2]);
    assert(3 == hist.nnonemptybins);
    assert(1.2 == hist.mins[0]);
    assert(2.3 == hist.mins[1]);
    assert(3.4 == hist.mins[2]);
    assert(2.3 == hist.maxs[0]);
    assert(3.4 == hist.maxs[1]);
    assert(4.5 == hist.maxs[2]);
    assert(0.45 == hist.percentinrange);
    assert(frequencies != hist.buffer);
    assert(1 == hist.buffer[0]);
    assert(9 == hist.buffer[1]);
    assert(4 == hist.buffer[2]);
    assert(2 == hist.buffer[3]);
    assert(6 == hist.buffer[4]);
    assert(50 == hist.buffer[5]);
}

void test_frequencies_to_histogram() {
    int32_t densevalues[] = {0, 0, 0, 4, 5, 6, 7, 8};
    int32_t nbins[] = {2, 2, 2};
    double mins[] = {0.1, 0.2, 0.3};
    double maxs[] = {0.4, 0.5, 0.6};
    Histogram densehist = frequencies_to_histogram(
            3, densevalues, nbins, mins, maxs, 0.98);
    assert(0 == densehist.issparse);
    int32_t sparsevalues[] = {0, 12, 56, 0, 0, 34, 2, 0};
    Histogram sparsehist = frequencies_to_histogram(
            3, sparsevalues, nbins, mins, maxs, 0.87);
    assert(1 == sparsehist.issparse);
}

void test_adjust_range_by_methods() {
    int32_t ngridx = 10, ngridy = 10, ngridz = 10;
    int32_t ngrids = ngridx * ngridy * ngridz;
    int32_t igrid, x, y, z;
    SamplingRegion samplingregion;
    char* methods[3];
    double vala = 0.0, valb = 0.0, valc = 990.0;
    double mins[3], maxs[3], adjmins[3], adjmaxs[3];
    samplingregion.xbeg = 3;
    samplingregion.xend = 8;
    samplingregion.ybeg = 2;
    samplingregion.yend = 4;
    samplingregion.zbeg = 0;
    samplingregion.zend = 10;
    samplingregion.ngridx = ngridx;
    samplingregion.ngridy = ngridy;
    samplingregion.ngridz = ngridz;
    samplingregion.ndims = 3;
    samplingregion.valuearrays[0] = malloc(ngrids * sizeof(double));
    samplingregion.valuearrays[1] = malloc(ngrids * sizeof(double));
    samplingregion.valuearrays[2] = malloc(ngrids * sizeof(double));
    for (z = samplingregion.zbeg; z < samplingregion.zend; ++z)
    for (y = samplingregion.ybeg; y < samplingregion.yend; ++y)
    for (x = samplingregion.xbeg; x < samplingregion.xend; ++x) {
        igrid = x + y * ngridx + z * ngridy * ngridx;
        samplingregion.valuearrays[0][igrid] = vala;
        samplingregion.valuearrays[1][igrid] = valb;
        samplingregion.valuearrays[2][igrid] = valc;
        vala += 0.1;
        valb += 1.0;
        valc -= 10.0;
    }
    methods[0] = malloc(30 * sizeof(char));
    strncpy(methods[0], "range", 5);
    methods[1] = malloc(30 * sizeof(char));
    strncpy(methods[1], "normalized_range", 16);
    methods[2] = malloc(30 * sizeof(char));
    strncpy(methods[2], "percent_range", 13);
    mins[0] = 0.5; maxs[0] = 0.7;
    mins[1] = 0.2; maxs[1] = 0.8;
    mins[2] = 0.1; maxs[2] = 0.9;

    adjust_range_by_methods(samplingregion.ndims, samplingregion, methods, mins,
            maxs, adjmins, adjmaxs);
    assert(adjmins[0] == 0.5 && adjmaxs[0] == 0.7);
    assert(fabs(adjmins[1] - 0.2 * 99.0) < 0.0001);
    assert(fabs(adjmaxs[1] - 0.8 * 99.0) < 0.0001);
    assert(fabs(adjmins[2] - 100.0) < 0.0001);
    assert(fabs(adjmaxs[2] - 900.0) < 0.0001);

    free(methods[0]);
    free(methods[1]);
    free(methods[2]);
    free(samplingregion.valuearrays[0]);
    free(samplingregion.valuearrays[1]);
    free(samplingregion.valuearrays[2]);
}

void test_compute_nbins() {
    char* nbinstrs[3];
    int32_t nbins[3];
    int32_t ngridx = 10, ngridy = 10, ngridz = 10;
    int32_t ngrids = ngridx * ngridy * ngridz;
    int32_t igrid, x, y, z;
    double vala = 0.0, valb = 0.0, valc = 990.0;
    double mins[3], maxs[3];
    SamplingRegion samplingregion;
    samplingregion.xbeg = 3;
    samplingregion.xend = 8;
    samplingregion.ybeg = 2;
    samplingregion.yend = 4;
    samplingregion.zbeg = 0;
    samplingregion.zend = 10;
    samplingregion.ngridx = ngridx;
    samplingregion.ngridy = ngridy;
    samplingregion.ngridz = ngridz;
    samplingregion.ndims = 3;
    samplingregion.valuearrays[0] = malloc(ngrids * sizeof(double));
    samplingregion.valuearrays[1] = malloc(ngrids * sizeof(double));
    samplingregion.valuearrays[2] = malloc(ngrids * sizeof(double));
    for (z = samplingregion.zbeg; z < samplingregion.zend; ++z)
    for (y = samplingregion.ybeg; y < samplingregion.yend; ++y)
    for (x = samplingregion.xbeg; x < samplingregion.xend; ++x) {
        igrid = x + y * ngridx + z * ngridy * ngridx;
        samplingregion.valuearrays[0][igrid] = vala;
        samplingregion.valuearrays[1][igrid] = valb;
        samplingregion.valuearrays[2][igrid] = valc;
        vala += 0.1;
        valb += 1.0;
        valc -= 10.0;
    }
    nbinstrs[0] = malloc(30 * sizeof(char));
    nbinstrs[1] = malloc(30 * sizeof(char));
    nbinstrs[2] = malloc(30 * sizeof(char));

    mins[0] = 0.999; maxs[0] = 2.0;
    mins[1] = 20.0; maxs[1] = 30.0;
    mins[2] = 200.0; maxs[2] = 300.0;
    strncpy(nbinstrs[0], "sturges", 7);
    strncpy(nbinstrs[1], "freedman", 8);
    strncpy(nbinstrs[2], "8\0", 2);
    compute_nbins(
            samplingregion.ndims, samplingregion, mins, maxs, nbinstrs, nbins);
    // printf("nbins[0] = %d, nbins[1] = %d, and nbins[2] = %d\n", nbins[0], nbins[1], nbins[2]);
    assert(nbins[0] == 5);
    assert(nbins[1] == 2);
    assert(nbins[2] == 8);

    mins[0] = 0.0; maxs[0] = 0.0;
    strncpy(nbinstrs[0], "sturges", 7);
    compute_nbins(1, samplingregion, mins, maxs, nbinstrs, nbins);
    assert(nbins[0] == 1);

    mins[0] = 5.0; maxs[0] = 0.0;
    strncpy(nbinstrs[0], "freedman", 8);
    compute_nbins(1, samplingregion, mins, maxs, nbinstrs, nbins);
    assert(nbins[0] == 1);

    mins[0] = 0.999; maxs[0] = 2.0;
    strncpy(nbinstrs[0], "incorrect", 9);
    compute_nbins(1, samplingregion, mins, maxs, nbinstrs, nbins);
    assert(nbins[0] == 2);

    free(nbinstrs[0]);
    free(nbinstrs[1]);
    free(nbinstrs[2]);
    free(samplingregion.valuearrays[0]);
    free(samplingregion.valuearrays[1]);
    free(samplingregion.valuearrays[2]);
}

void test_generate_histogram_from_sampling_region() {
    int32_t i, nbins[1];
    char* methods[1];
    double mins[1], maxs[1];
    SamplingRegion samplingregion;
    samplingregion.xbeg = 3;
    samplingregion.xend = 8;
    samplingregion.ybeg = 2;
    samplingregion.yend = 4;
    samplingregion.zbeg = 0;
    samplingregion.zend = 9;
    samplingregion.ngridx = 10;
    samplingregion.ngridy = 10;
    samplingregion.ngridz = 10;
    samplingregion.ndims = 1;
    samplingregion.valuearrays[0] = malloc(
            samplingregion.ngridx * samplingregion.ngridy *
            samplingregion.ngridz * sizeof(double));
    for (i = 0; i < samplingregion.ngridx * samplingregion.ngridy * samplingregion.ngridz; ++i) {
        samplingregion.valuearrays[0][i] = 5.0;
    }
    nbins[0] = 20;
    mins[0] = 4.0;
    maxs[0] = 7.0;
    methods[0] = malloc(30);
    strncpy(methods[0], "range", 5);
    Histogram hist = generate_histogram_from_sampling_region(
            1, samplingregion, nbins, mins, maxs);
    free(methods[0]);
    assert(1 == hist.ndims);
    assert(1 == hist.nnonemptybins);
    assert(6 == hist.buffer[0]);
    assert(90 == hist.buffer[1]);
}

void test_write_histogram_meta() {
    FILE* file, *read;
    int32_t ndims, ngridx, ngridy, ngridz, nhistx, nhisty, nhistz;
    int32_t nbins[] = {5, 6, 7};
    double log_base[] = {0.2, 0.3, 0.4};
    file = fopen("testfile.bin", "wb");
    write_histogram_meta(file, 4, 5, 6, 7, 8, 9, 3, nbins, log_base);
    fclose(file);
    read = fopen("testfile.bin", "rb");
    fread(&ndims, 1, sizeof(int32_t), read);
    fread(&ngridx, 1, sizeof(int32_t), read);
    fread(&ngridy, 1, sizeof(int32_t), read);
    fread(&ngridz, 1, sizeof(int32_t), read);
    fread(&nhistx, 1, sizeof(int32_t), read);
    fread(&nhisty, 1, sizeof(int32_t), read);
    fread(&nhistz, 1, sizeof(int32_t), read);
    fread(log_base, ndims, sizeof(double), read);
    fclose(read);
    assert(3 == ndims);
    assert(4 == ngridx);
    assert(5 == ngridy);
    assert(6 == ngridz);
    assert(7 == nhistx);
    assert(8 == nhisty);
    assert(9 == nhistz);
    assert(fabs(0.2 - log_base[0]) < 0.0001);
    assert(fabs(0.3 - log_base[1]) < 0.0001);
    assert(fabs(0.4 - log_base[2]) < 0.0001);
}

void test_write_histogram_1d() {
    FILE *write, *read;
    Histogram writehist, readhist;
    int32_t i;
    int32_t writebuffer[5] = {0, 0, 1, 2, 3};
    int32_t readbuffer[5];
    writehist.ndims = 1;
    writehist.buffer = writebuffer;
    writehist.issparse = 0;
    writehist.nbins[0] = 5;
    writehist.nnonemptybins = 3;
    writehist.mins[0] = 0.2;
    writehist.maxs[0] = 2.3;
    writehist.percentinrange = 0.63;
    writehist.nnonemptybins = 3;
    write = fopen("testfile.bin", "wb");
    write_histogram(write, writehist);
    fclose(write);
    read = fopen("testfile.bin", "rb");
    readhist.buffer = readbuffer;
    fread(&readhist.issparse, 1, sizeof(int32_t), read);
    fread(&readhist.mins[0], 1, sizeof(double), read);
    fread(&readhist.maxs[0], 1, sizeof(double), read);
    fread(&readhist.nbins[0], 1, sizeof(int32_t), read);
    fread(&readhist.percentinrange, 1, sizeof(double), read);
    fread(&readhist.nnonemptybins, 1, sizeof(int32_t), read);
    fread(readhist.buffer, 5, sizeof(int32_t), read);
    fclose(read);
    assert(0 == readhist.issparse);
    assert(0.2 == readhist.mins[0]);
    assert(2.3 == readhist.maxs[0]);
    assert(5 == readhist.nbins[0]);
    assert(0.63 == readhist.percentinrange);
    assert(3 == readhist.nnonemptybins);
    for (i = 0; i < 5; ++i)
        assert(writebuffer[i] == readbuffer[i]);
}

void test_write_histogram_3d() {
    FILE *write, *read;
    Histogram writehist, readhist;
    int32_t i;
    int32_t writebuffer[8] = {0, 0, 1, 2, 3, 0, 5, 6};
    int32_t readbuffer[8];
    writehist.ndims = 3;
    writehist.buffer = writebuffer;
    writehist.issparse = 0;
    writehist.nbins[0] = 2;
    writehist.nbins[1] = 2;
    writehist.nbins[2] = 2;
    writehist.nnonemptybins = 3;
    writehist.mins[0] = 0.2;
    writehist.mins[1] = 0.3;
    writehist.mins[2] = 0.4;
    writehist.maxs[0] = 2.3;
    writehist.maxs[1] = 2.4;
    writehist.maxs[2] = 2.5;
    writehist.percentinrange = 0.63;
    writehist.nnonemptybins = 3;
    write = fopen("testfile.bin", "wb");
    write_histogram(write, writehist);
    fclose(write);
    read = fopen("testfile.bin", "rb");
    readhist.buffer = readbuffer;
    fread(&readhist.issparse, 1, sizeof(int32_t), read);
    fread(&readhist.mins[0], 3, sizeof(double), read);
    fread(&readhist.maxs[0], 3, sizeof(double), read);
    fread(&readhist.nbins[0], 3, sizeof(int32_t), read);
    fread(&readhist.percentinrange, 1, sizeof(double), read);
    fread(&readhist.nnonemptybins, 1, sizeof(int32_t), read);
    fread(readhist.buffer, 8, sizeof(int32_t), read);
    fclose(read);
    assert(0 == readhist.issparse);
    assert(0.2 == readhist.mins[0]);
    assert(0.3 == readhist.mins[1]);
    assert(0.4 == readhist.mins[2]);
    assert(2.3 == readhist.maxs[0]);
    assert(2.4 == readhist.maxs[1]);
    assert(2.5 == readhist.maxs[2]);
    assert(2 == readhist.nbins[0]);
    assert(2 == readhist.nbins[1]);
    assert(2 == readhist.nbins[2]);
    assert(0.63 == readhist.percentinrange);
    assert(3 == readhist.nnonemptybins);
    for (i = 0; i < 8; ++i)
        assert(writebuffer[i] == readbuffer[i]);
}

void test_serialize_domain_histograms() {
    double log_base[] = {0.2};
    char *buffer;
    int32_t nbyte, ibyte, ngridx, ngridy, ngridz, nhistx, nhisty, nhistz, ndims;
    double read_log_base[1];
    Histogram writehists[2];
    Histogram readhist1, readhist2;
    int32_t writebuffer1[3] = {3, 4, 5}, writebuffer2[4] = {6, 7, 8, 9};
    int32_t readbuffer1[3], readbuffer2[4];
    // write histogram 1
    writehists[0].ndims = 1;
    writehists[0].buffer = writebuffer1;
    writehists[0].issparse = 0;
    writehists[0].nbins[0] = 3;
    writehists[0].nnonemptybins = 3;
    writehists[0].mins[0] = 0.2;
    writehists[0].mins[1] = 0.3;
    writehists[0].mins[2] = 0.4;
    writehists[0].maxs[0] = 2.2;
    writehists[0].maxs[1] = 2.3;
    writehists[0].maxs[2] = 2.4;
    writehists[0].percentinrange = 0.63;
    // write histogram 2
    writehists[1].ndims = 1;
    writehists[1].buffer = writebuffer2;
    writehists[1].issparse = 0;
    writehists[1].nbins[0] = 4;
    writehists[1].nnonemptybins = 4;
    writehists[1].mins[0] = 0.5;
    writehists[1].mins[1] = 0.6;
    writehists[1].mins[2] = 0.7;
    writehists[1].maxs[0] = 2.5;
    writehists[1].maxs[1] = 2.6;
    writehists[1].maxs[2] = 2.7;
    writehists[1].percentinrange = 0.96;
    serialize_domain_histograms(4, 2, 2, 2, 1, 1, 1, log_base, writehists, 2,
            &buffer, &nbyte);
    ibyte = 0;
    // read domain meta
    memcpy(&ndims, buffer + ibyte, sizeof(int32_t)); ibyte += sizeof(int32_t);
    memcpy(&ngridx, buffer + ibyte, sizeof(int32_t)); ibyte += sizeof(int32_t);
    memcpy(&ngridy, buffer + ibyte, sizeof(int32_t)); ibyte += sizeof(int32_t);
    memcpy(&ngridz, buffer + ibyte, sizeof(int32_t)); ibyte += sizeof(int32_t);
    memcpy(&nhistx, buffer + ibyte, sizeof(int32_t)); ibyte += sizeof(int32_t);
    memcpy(&nhisty, buffer + ibyte, sizeof(int32_t)); ibyte += sizeof(int32_t);
    memcpy(&nhistz, buffer + ibyte, sizeof(int32_t)); ibyte += sizeof(int32_t);
    memcpy(read_log_base, buffer + ibyte, sizeof(double) * ndims); ibyte += sizeof(double) * ndims;
    // read histogram 1
    readhist1.buffer = readbuffer1;
    readhist1.ndims = ndims;
    memcpy(&readhist1.issparse, buffer + ibyte, sizeof(int32_t)); ibyte += sizeof(int32_t);
    memcpy(readhist1.mins, buffer + ibyte, sizeof(double) * ndims); ibyte += sizeof(double) * ndims;
    memcpy(readhist1.maxs, buffer + ibyte, sizeof(double) * ndims); ibyte += sizeof(double) * ndims;
    memcpy(readhist1.nbins, buffer + ibyte, sizeof(int32_t) * ndims); ibyte += sizeof(int32_t) * ndims;
    memcpy(&readhist1.percentinrange, buffer + ibyte, sizeof(double)); ibyte += sizeof(double);
    memcpy(&readhist1.nnonemptybins, buffer + ibyte, sizeof(int32_t)); ibyte += sizeof(int32_t);
    memcpy(readhist1.buffer, buffer + ibyte, sizeof(int32_t) * readhist1.nbins[0]); ibyte += sizeof(int32_t) * readhist1.nbins[0];
    // read histogram 2
    readhist2.buffer = readbuffer2;
    readhist2.ndims = ndims;
    memcpy(&readhist2.issparse, buffer + ibyte, sizeof(int32_t)); ibyte += sizeof(int32_t);
    memcpy(readhist2.mins, buffer + ibyte, sizeof(double) * ndims); ibyte += sizeof(double) * ndims;
    memcpy(readhist2.maxs, buffer + ibyte, sizeof(double) * ndims); ibyte += sizeof(double) * ndims;
    memcpy(readhist2.nbins, buffer + ibyte, sizeof(int32_t) * ndims); ibyte += sizeof(int32_t) * ndims;
    memcpy(&readhist2.percentinrange, buffer + ibyte, sizeof(double)); ibyte += sizeof(double);
    memcpy(&readhist2.nnonemptybins, buffer + ibyte, sizeof(int32_t)); ibyte += sizeof(int32_t);
    memcpy(readhist2.buffer, buffer + ibyte, sizeof(int32_t) * readhist2.nbins[0]); ibyte += sizeof(int32_t) * readhist2.nbins[0];
    // assert meta
    assert(ndims == 1);
    assert(ngridx == 4);
    assert(ngridy == 2);
    assert(ngridz == 2);
    assert(nhistx == 2);
    assert(nhisty == 1);
    assert(nhistz == 1);
    assert(log_base[0] == read_log_base[0]);
    // assert histogram 1
    assert(writehists[0].ndims == readhist1.ndims);
    assert(writehists[0].issparse == readhist1.issparse);
    assert(writehists[0].mins[0] == readhist1.mins[0]);
    assert(writehists[0].maxs[0] == readhist1.maxs[0]);
    assert(writehists[0].nbins[0] == readhist1.nbins[0]);
    assert(writehists[0].percentinrange == readhist1.percentinrange);
    assert(writehists[0].nnonemptybins == readhist1.nnonemptybins);
    assert(writehists[0].buffer[0] == readhist1.buffer[0]);
    assert(writehists[0].buffer[1] == readhist1.buffer[1]);
    assert(writehists[0].buffer[2] == readhist1.buffer[2]);
    // assert histogram 2
    assert(writehists[1].ndims == readhist2.ndims);
    assert(writehists[1].issparse == readhist2.issparse);
    assert(writehists[1].mins[0] == readhist2.mins[0]);
    assert(writehists[1].maxs[0] == readhist2.maxs[0]);
    assert(writehists[1].nbins[0] == readhist2.nbins[0]);
    assert(writehists[1].percentinrange == readhist2.percentinrange);
    assert(writehists[1].nnonemptybins == readhist2.nnonemptybins);
    assert(writehists[1].buffer[0] == readhist2.buffer[0]);
    assert(writehists[1].buffer[1] == readhist2.buffer[1]);
    assert(writehists[1].buffer[2] == readhist2.buffer[2]);
    assert(writehists[1].buffer[3] == readhist2.buffer[3]);
    free(buffer);
}

// void test_write_domain_histograms() {
//     FILE* file, *read;
//     int32_t ndims, ngridx, ngridy, ngridz, nhistx, nhisty, nhistz;
//     int32_t nbins[] = {5, 6, 7};
//     double log_base[] = {0.2, 0.3, 0.4};
//     file = fopen("testfile.bin", "wb");

//     write_histogram_meta(file, 4, 5, 6, 7, 8, 9, 3, nbins, log_base);
//     fclose(file);
//     read = fopen("testfile.bin", "rb");
//     fread(&ndims, 1, sizeof(int32_t), read);
//     fread(&ngridx, 1, sizeof(int32_t), read);
//     fread(&ngridy, 1, sizeof(int32_t), read);
//     fread(&ngridz, 1, sizeof(int32_t), read);
//     fread(&nhistx, 1, sizeof(int32_t), read);
//     fread(&nhisty, 1, sizeof(int32_t), read);
//     fread(&nhistz, 1, sizeof(int32_t), read);
//     fread(log_base, ndims, sizeof(double), read);
//     fclose(read);
//     assert(3 == ndims);
//     assert(4 == ngridx);
//     assert(5 == ngridy);
//     assert(6 == ngridz);
//     assert(7 == nhistx);
//     assert(8 == nhisty);
//     assert(9 == nhistz);
//     assert(fabs(0.2 - log_base[0]) < 0.0001);
//     assert(fabs(0.3 - log_base[1]) < 0.0001);
//     assert(fabs(0.4 - log_base[2]) < 0.0001);
// }

void test_ids_to_flat() {
    int32_t dims[] = {5, 6, 7};
    int32_t ids[] = {2, 3, 4};
    assert(2 + 15 + 120 == ids_to_flat(3, dims, ids));
}

void test_values_to_bin_index() {
    double values[] = {2.0, 5.0, 7.5};
    double mins[] = {0.0, 0.0, 0.0};
    double maxs[] = {10.0, 10.0, 10.0};
    int32_t nbins[] = {10, 10, 10};
    assert(2 + 50 + 700 == values_to_bin_index(3, values, mins, maxs, nbins));
}

int main(void) {
    test_total_number_of_bins();
    test_get_number_of_nonempty_bins();
    test_frequencies_to_histogram_1d_dense();
    test_frequencies_to_histogram_2d_dense();
    test_frequencies_to_histogram_3d_dense();
    test_frequencies_to_histogram_1d_sparse();
    test_frequencies_to_histogram_2d_sparse();
    test_frequencies_to_histogram_3d_sparse();
    test_frequencies_to_histogram();
    test_adjust_range_by_methods();
    test_compute_nbins();
    test_generate_histogram_from_sampling_region();
    test_write_histogram_meta();
    test_write_histogram_1d();
    test_write_histogram_3d();
    test_serialize_domain_histograms();
    test_ids_to_flat();
    test_values_to_bin_index();
    return 0;
}
