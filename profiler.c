#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <mpi.h>
#include "profiler.h"

//
//
//
// internal variables
//
//
//

// constant variables
static const char filename[] = "../data/timing.csv";
#define maxtitles 100
#define maxrows 20
// the data which will be written to file
static char titles[maxtitles][100];
static double times[maxtitles][maxrows];
static double timesteps[maxrows];
// the stopwatches
static double timestarts[maxtitles];
// variables to keep track of the state
static int ntitles = 0;
static int nrows = 0;
static int currrow = 0;

//
//
//
// internal functions
//
//
//

static void cleartitles() {
    // temp variables
    int ititle = 0;
    // clear the data buffers
    for (ititle = 0; ititle < maxtitles; ++ititle) {
        titles[ititle][0] = '\0';
    }
    // reset the state variables
    ntitles = 0;
}

static void clearbuffer() {
    // temp variables
    int ititle = 0;
    int irow = 0;
    // clear the data buffers
    for (ititle = 0; ititle < maxtitles; ++ititle)
    for (irow = 0; irow < maxrows; ++irow) {
        times[ititle][irow] = 0.0;
    }
    for (irow = 0; irow < maxrows; ++irow) {
        timesteps[irow] = 0.0;
    }
    // reset the state variables
    nrows = 0;
    currrow = 0;
}

static void clearfile() {
    fclose(fopen(filename, "w"));
}

static int t2i(const char* title) {
    int ititle = 0;
    for (ititle = 0; ititle < ntitles; ++ititle) {
        if (0 == strcmp(title, titles[ititle])) {
            return ititle;
        }
    }
    // if exceeds the maximum allowed titles
    if (ititle >= maxtitles)
        return -1;
    // if title doesn't exist yet, create it
    strcpy(titles[ititle], title);
    ++ntitles;
    return ititle;
}

static double now() {
    double ret;
    struct timeval t;
    gettimeofday(&t, NULL);
    ret = t.tv_sec + t.tv_usec / 1000000.0;
    return ret;
}

void reset(int ititle) {
    times[ititle][currrow] = 0.0;
}

void start(int ititle) {
    timestarts[ititle] = now();
}

void pause(int ititle) {
    // printf("pause: %d, %d, %f, %f\n", ititle, currrow, times[ititle][currrow], now()-timestarts[ititle]);
    times[ititle][currrow] += (now() - timestarts[ititle]);
}

void writetitles(FILE* pfile) {
    int ititle = 0;
    if (ntitles <= 0)
        return;
    fprintf(pfile, "time");
    for (ititle = 0; ititle < ntitles; ++ititle) {
        fprintf(pfile, ",%s", titles[ititle]);
    }
    fprintf(pfile, "\n");
}

void writebuffer(FILE* pfile) {
    int ititle = 0;
    int irow = 0;
    for (irow = 0; irow < nrows; ++irow) {
        fprintf(pfile, "%E", timesteps[irow]);
        for (ititle = 0; ititle < ntitles; ++ititle) {
            // printf("ititle: %d, irow: %d, time: %f\n", ititle, irow, times[ititle][irow]);
            fprintf(pfile, ",%f", times[ititle][irow]);
        }
        fprintf(pfile, "\n");
    }
}

//
//
//
// the public interface implementation
//
//
//

void profiler_clear_() {
    cleartitles();
    clearbuffer();
    clearfile();
}

// TODO: refactor using function pointer
void profiler_advance_(const double* ptimestep) {
    int ititle = 0;
    int nextrow = currrow + 1;
    // set time step
    timesteps[currrow] = *ptimestep;
    // output to file if the buffer is full
    if (nextrow >= maxrows) {
        profiler_appendfile_();
        clearbuffer();
        nextrow = 0;
    }
    // update current row variable
    currrow = nextrow;
    nrows = currrow;
    // and reset the stopwatch
    timestarts[currrow] = 0.0;
    for (ititle = 0; ititle < maxtitles; ++ititle) {
        reset(ititle);
    }
}

void profiler_mpi_advance_(const double* ptimestep, int* pmyid, int* proot, MPI_Fint* pcomm) {
    int ititle = 0;
    int nextrow = currrow + 1;
    // set time step
    timesteps[currrow] = *ptimestep;
    // output to file if the buffer is full
    if (nextrow >= maxrows) {
        profiler_mpi_appendfile_(pmyid, proot, pcomm);
        clearbuffer();
        nextrow = 0;
    }
    // update current row variable
    currrow = nextrow;
    nrows = currrow;
    // and reset the stopwatch
    timestarts[currrow] = 0.0;
    for (ititle = 0; ititle < maxtitles; ++ititle) {
        reset(ititle);
    }
}

void profiler_reset_(const char* title) {
    int ititle = t2i(title);
    reset(ititle);
}

void profiler_start_(const char* title) {
    int ititle = t2i(title);
    start(ititle);
}

void profiler_pause_(const char* title) {
    int ititle = t2i(title);
    pause(ititle);
}

void profiler_appendfile_() {
    long size = 0;
    // open file
    FILE* pfile = fopen(filename, "a");
    // to the end of the file
    fseek(pfile, 0, SEEK_END);
    // if file is empty, write the titles
    size = ftell(pfile);
    if (0 == size) {
        writetitles(pfile);
    }
    // append the buffer to the end of the file
    writebuffer(pfile);
    // close the file
    fclose(pfile);   
}

void profiler_mpi_appendfile_(int* pmyid, int* proot, MPI_Fint* pcomm) {
    int myid = *pmyid;
    int root = *proot;
    MPI_Comm comm = MPI_Comm_f2c(*pcomm);
    double recv[maxtitles][maxrows];
    // reduce to the root process so only 1 file is being output
    MPI_Reduce((void*)times, (void*)recv,
               maxtitles * maxrows, MPI_DOUBLE,
               MPI_MAX, root, comm);
    memcpy((void*)times, (void*)recv, maxtitles * maxrows * sizeof(double));
    // output
    if (myid == root) {
        profiler_appendfile_();
    }
}
