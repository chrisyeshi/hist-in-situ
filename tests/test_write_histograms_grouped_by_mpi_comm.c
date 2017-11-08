#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/stat.h>
#include <mpi.h>
#include <histogram.h>
#include <histogram_in.h>

#define nDims 1
#define iConfig 1
#define rootId 0
#define timestep 0.001

void create_mpi_comms(int nComms, MPI_Comm* comm, int* groupId) {
    int world_rank, world_size, color;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    *groupId = world_rank % nComms;
    MPI_Comm_split(MPI_COMM_WORLD, *groupId, world_rank, comm);
}

Histogram random_histogram() {
        // int32_t frequencies[5] = {0, 0, 1, 2, 3};
        // int32_t nbins[] = {5}, nBytes, histBufferByteCount;
        // double mins[] = {12.34};
        // double maxs[] = {34.56};
        // Histogram hist = frequencies_to_histogram_dense(
        //         1, frequencies, nbins, 3, mins, maxs, 0.78);


    int iBin;
    int32_t frequencies[5] = {0, 0, 1, 2, 3};
    int32_t nbins[] = {5};
    double mins[] = {12.34};
    double maxs[] = {34.56};
    // for (int iBin = 0; iBin < 5; ++iBin) {
    //     frequencies[iBin] = rand();
    // }
    Histogram hist = frequencies_to_histogram_dense(
            nDims, frequencies, nbins, 3, mins, maxs, 0.78);
    return hist;
}

void create_histograms(int nHistsPerRank, Histogram hists[]) {
    int iHist;
    for (iHist = 0; iHist < nHistsPerRank; ++iHist) {
        hists[iHist] = random_histogram();
    }
}

void write_histograms(
        MPI_Comm comm, int nHistsPerRank, Histogram hists[], int groupId) {
    double log_base[] = {0.0};
    int32_t rankInComm, commSize;
    MPI_Comm_rank(comm, &rankInComm);
    MPI_Comm_size(comm, &commSize);
    mkdir("../data_pdf/", S_IRWXU);
    mkdir("../data_pdf/pdf-1.0000E-03/", S_IRWXU);
    write_histograms_grouped_by_mpi_comm(
            5 /*ngridx*/, 5 /*ngridy*/, 5 /*ngridz*/,
            1 /*nhistx*/, 1 /*nhisty*/, nHistsPerRank /*nhistz*/,
            nDims, log_base, hists, nHistsPerRank /*nhists*/,
            timestep, rankInComm /*yid*/, 0 /*rootid*/,
            commSize /*ypes*/, comm, groupId, iConfig, "test" /*midname*/);
}

int is_root(MPI_Comm comm) {
    int rank;
    MPI_Comm_rank(comm, &rank);
    return rank == rootId;
}

void gather_histograms(int nHistsPerRank, Histogram hists[],
        Histogram** outHists, int* nHists) {
    char *bufferPerRank, *buffer;
    int nBytesPerRank, nBytes, *nBytesPerRanks;

    serialize_histograms(nHistsPerRank, hists, &nBytesPerRank, &bufferPerRank);
    gather_buffers(nBytesPerRank, bufferPerRank, &nBytesPerRanks, &nBytes,
            &buffer, rootId, MPI_COMM_WORLD);
    if (is_root(MPI_COMM_WORLD)) {
        deserialize_histograms(nBytes, buffer, nDims, nHists, outHists);
    }

    free(bufferPerRank);
    if (is_root(MPI_COMM_WORLD)) {
        free(nBytesPerRanks);
        free(buffer);
    }
}

void read_histograms(Histogram** hists, int* nHists) {
    FILE *f1, *f2;
    char *buffer;
    long fsize1, fsize2;
    int nBytes = 560, offset = 0;
    int metaByteCount = 36, histByteCountPerRank = 112;
    int iHist;
    f1 = fopen("../data_pdf/pdf-1.0000E-03/pdfs-test-001.00000", "rb");
    f2 = fopen("../data_pdf/pdf-1.0000E-03/pdfs-test-001.00001", "rb");
    buffer = malloc(nBytes);

    fseek(f1, metaByteCount, SEEK_CUR);
    fread(buffer + offset, histByteCountPerRank, 1, f1);
    offset += histByteCountPerRank;

    fseek(f1, metaByteCount, SEEK_CUR);
    fread(buffer + offset, histByteCountPerRank, 1, f1);
    offset += histByteCountPerRank;

    fseek(f1, metaByteCount, SEEK_CUR);
    fread(buffer + offset, histByteCountPerRank, 1, f1);
    offset += histByteCountPerRank;

    fseek(f2, metaByteCount, SEEK_CUR);
    fread(buffer + offset, histByteCountPerRank, 1, f2);
    offset += histByteCountPerRank;

    fseek(f2, metaByteCount, SEEK_CUR);
    fread(buffer + offset, histByteCountPerRank, 1, f2);
    offset += histByteCountPerRank;

    fclose(f1);
    fclose(f2);

    deserialize_histograms(nBytes, buffer, nDims, nHists, hists);
}

void compare_histograms(int inHistCount, Histogram inHists[], int outHistCount,
        Histogram outHists[]) {
    printf("%d : %d\n", inHistCount, outHistCount);
    assert(inHistCount == outHistCount);
}

void deallocate_histograms(int nHists, Histogram hists[]) {
    int iHist;
    for (iHist = 0; iHist < nHists; ++iHist) {
        deallocate_histogram(hists[iHist]);
    }
}

void test_histograms(int nHistsPerRank, Histogram hists[]) {
    Histogram *outHists, *inHists;
    int outHistCount, inHistCount;
    gather_histograms(nHistsPerRank, hists, &outHists, &outHistCount);
    if (is_root(MPI_COMM_WORLD)) {
        read_histograms(&inHists, &inHistCount);
        compare_histograms(inHistCount, inHists, outHistCount, outHists);
        deallocate_histograms(inHistCount, inHists);
        free(inHists);
        deallocate_histograms(outHistCount, outHists);
        free(outHists);
    }
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    int groupId, worldRank;
    MPI_Comm comm;
    Histogram hists[2];

    create_mpi_comms(2, &comm, &groupId);
    create_histograms(2, hists);
    write_histograms(comm, 2, hists, groupId);
    test_histograms(2, hists);
    deallocate_histograms(2, hists);

    MPI_Comm_free(&comm);

    MPI_Finalize();
    return 0;
}