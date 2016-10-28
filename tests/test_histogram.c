#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <histogram.h>

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

void test_frequencies_to_histogram_dense() {
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

void test_frequencies_to_histogram_sparse() {
    int32_t frequencies[4] = {0, 5, 0, 0};
    int32_t nbins[] = {2, 2};
    double mins[] = {12.34};
    double maxs[] = {34.56};
    Histogram hist = frequencies_to_histogram_sparse(
            2, frequencies, nbins, 1, mins, maxs, 0.78);
    assert(2 == hist.ndims);
    assert(1 == hist.issparse);
    assert(2 == hist.nbins[0]);
    assert(2 == hist.nbins[1]);
    assert(1 == hist.nnonemptybins);
    assert(12.34 == hist.mins[0]);
    assert(34.56 == hist.maxs[0]);
    assert(0.78 == hist.percentinrange);
    assert(frequencies != hist.buffer);
    assert(1 == hist.buffer[0]);
    assert(5 == hist.buffer[1]);
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

void test_generate_histogram_from_sampling_region() {
    int32_t i, nbins[1];
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
    Histogram hist = generate_histogram_from_sampling_region(
            1, samplingregion, nbins, mins, maxs);
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
    fread(nbins, ndims, sizeof(int32_t), read);
    fread(log_base, ndims, sizeof(double), read);
    fclose(read);
    assert(3 == ndims);
    assert(4 == ngridx);
    assert(5 == ngridy);
    assert(6 == ngridz);
    assert(7 == nhistx);
    assert(8 == nhisty);
    assert(9 == nhistz);
    assert(5 == nbins[0]);
    assert(6 == nbins[1]);
    assert(7 == nbins[2]);
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
    fread(&readhist.percentinrange, 1, sizeof(double), read);
    fread(&readhist.nnonemptybins, 1, sizeof(int32_t), read);
    fread(readhist.buffer, 5, sizeof(int32_t), read);
    fclose(read);
    assert(0 == readhist.issparse);
    assert(0.2 == readhist.mins[0]);
    assert(2.3 == readhist.maxs[0]);
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
    assert(0.63 == readhist.percentinrange);
    assert(3 == readhist.nnonemptybins);
    for (i = 0; i < 8; ++i)
        assert(writebuffer[i] == readbuffer[i]);
}

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
    test_frequencies_to_histogram_dense();
    test_frequencies_to_histogram_sparse();
    test_frequencies_to_histogram();
    test_generate_histogram_from_sampling_region();
    test_write_histogram_meta();
    test_write_histogram_1d();
    test_write_histogram_3d();
    test_ids_to_flat();
    test_values_to_bin_index();

    return 0;
}