#ifndef __PROFILER_H__
#define __PROFILER_H__

#include <mpi.h>

//
//
//
// the public interface
//
//
//

// clears both the titles and times arrays and the output file
void profiler_clear_();
// advance to next row in the time table, append to file if it's full
// the timestep is for indicating which time step the next row belongs to
void profiler_advance_(const double* ptimestep);
// mpi version of the advance function
void profiler_mpi_advance_(const double* ptimestep, int* pmyid, int* proot, MPI_Fint* pcomm);
// reset the stopwatch reading of the current time step
void profiler_reset_(const char* title);
// start the stopwatch of the current time step
void profiler_start_(const char* title);
// pause/stop the stopwatch of the current time step
void profiler_pause_(const char* title);
// finish timing and append to file
void profiler_appendfile_();
// mpi version of the appendfile function
void profiler_mpi_appendfile_(int* pmyid, int* proot, MPI_Fint* pcomm);

#endif // __PROFILER_H__