#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include "profiler.h"

//REMEMBER:
//The order of multidimensional array indexes is the opposite in fortran!
//Fortran starts indexing at 1!
//Make sure number of dom. decomp. and number of PDFs always divides evenly!

//TRACERS SORTED OUTPUT:
//organizes particles according to distribution function grid and saves particle data
void tracerssortedoutput_(double *ts, //the current timestep of the simulation
                          int32_t *myid, //the id of the parallel processing node
                          int32_t *nhx,
                          int32_t *nhy,
                          int32_t *nhz,
                          int32_t *nx, //number of cells in x direction
                          int32_t *ny, //number of cells in y direction
                          int32_t *nz, //number of cells in z direction
                          int64_t *ssn, //pointer to particle id array - ssn[fill]
                          double *iloc, //2d array of particle grid location - loc[3][ARRAY_SIZE]
                          int32_t *iARRAY_SIZE, //for indexing iloc
                          double *ixloc, //2d array of particle global location - xloc[3][fill]
                          double *isfield, //2d array of additional variable information - sfield[nfield][fill]
                          int32_t *fill, //number of particles
                          int32_t *nfield) //number of additional variables
{
    int i,j,k; //for loops
    FILE *write_ptr; //file output
    FILE *config_ptr;

    //delcare all arrays
    int *helper, *PDFindex, *sorter, *offsets;
    int64_t *ssn_out;
    double *xloc_out;
    double *sfield_out;
    char idname[100], helpername[100], xlocname[100], tracername[100];
    int ARRAY_SIZE = *iARRAY_SIZE;
    double (*loc)[ARRAY_SIZE] = (double (*)[ARRAY_SIZE]) iloc;
    double (*xloc)[*fill] = (double (*)[*fill]) ixloc;
    double (*sfield)[*fill] = (double (*)[*fill]) isfield;
    int32_t nPDFx, nPDFy, nPDFz, sizePDFx, sizePDFy, sizePDFz;

    nPDFx = (*nhx);
    nPDFy = (*nhy);
    nPDFz = (*nhz);

    // geometry of PDF regions (number in each direction and their size)
    sizePDFx = (*nx)/nPDFx;
    sizePDFy = (*ny)/nPDFy;
    sizePDFz = (*nz)/nPDFz;

    // begin profiler
    profiler_start_("sort particles");

    // helper file which contains start read locations and counts
    helper = (int *)calloc(2 * (nPDFx*nPDFy*nPDFz), sizeof(int)); //{start1, count1, start2, count2, etc...}
    
    //determine local row major PDF region index from particle position and store in array
    PDFindex = (int *)malloc(sizeof(int) * (*fill));    
    for(j=0; j < (*fill); j++) {
        int PDFPosx = (int)((((int)(loc[0][j])-1)%(*nx))/sizePDFx); //local PDF positions
        int PDFPosy = (int)((((int)(loc[1][j])-1)%(*ny))/sizePDFy);
        int PDFPosz = (int)((((int)(loc[2][j])-1)%(*nz))/sizePDFz);
        int index = PDFPosx + PDFPosy*nPDFx + PDFPosz*nPDFx*nPDFy;
        PDFindex[j] = index;
        helper[2*index + 1]++; //increment counts helper[i*2 + 1]
    }

    // determine start read locations helper[i*2 + 0]
    helper[0*2 + 0] = 0;
    for(i=1; i < nPDFx*nPDFy*nPDFz; ++i) {
        helper[i*2 + 0] = helper[(i-1)*2 + 0] + helper[(i-1)*2 + 1];
    }
    
    // build list of indexes of where each particle will go in final output buffers
    sorter = (int *)malloc(sizeof(int) * (*fill)); //stores the index of where each particle will go    
    offsets = (int *)calloc(nPDFx*nPDFy*nPDFz, sizeof(int)); //offsets for when we already placed a particle
    for(j=0; j < (*fill); j++) {
        int index = PDFindex[j]; //current pdf index
        sorter[j] = helper[index*2 + 0] + offsets[index]; //mark that particle j will move to start + offset
        offsets[index]++;
    }

    // output buffers
    ssn_out = (int64_t *)malloc(sizeof(int64_t) * (*fill)); //id (just one variable)
    xloc_out = (double *)malloc(sizeof(double) * 3 * (*fill)); //x,y,z (3 variables interleaved)
    sfield_out = (double *)malloc(sizeof(double) * (*nfield) * (*fill)); //all other variables (nfield variables NOT interleaved)

    for(j=0; j < (*fill); j++) {
        ssn_out[sorter[j]] = ssn[j]; //id's of particles
        xloc_out[3 * sorter[j] + 0] = xloc[0][j]; //physical particle positions
        xloc_out[3 * sorter[j] + 1] = xloc[1][j];
        xloc_out[3 * sorter[j] + 2] = xloc[2][j];
        for(k=0; k < (*nfield); k++) {
            sfield_out[k * (*fill) + sorter[j]] = sfield[k][j]; // other variables
        }
    }

    // switch profiler
    profiler_pause_("sort particles");
    profiler_start_("new output particles");

    // DomainSortedFormat
    sprintf(tracername, "../data/tracer-%.4E/pdfsortedtracer.%05d", *ts, *myid);
    write_ptr = fopen(tracername, "wb");
    // write the number of PDFs in each direction
    fwrite(&nPDFx, sizeof(int32_t), 1, write_ptr);
    fwrite(&nPDFy, sizeof(int32_t), 1, write_ptr);
    fwrite(&nPDFz, sizeof(int32_t), 1, write_ptr);
    // write the offsets and counts of particles for each PDF
    fwrite(helper, sizeof(int32_t), 2 * nPDFx * nPDFy * nPDFz, write_ptr);
    // write the particle ids
    fwrite(ssn_out, sizeof(int64_t), *fill, write_ptr);
    // write the particle physical positions
    fwrite(xloc_out, sizeof(double), 3 * (*fill), write_ptr);
    // write the particle scalar value
    fwrite(sfield_out, sizeof(double), (*nfield) * (*fill), write_ptr);
    fclose(write_ptr);




/*
    // ManyFilesFormat
    //write helper file
    sprintf(helpername, "../data/tracer-%.4E/sortedhelper.%05d", *ts, *myid);
    write_ptr = fopen(helpername, "wb");
    fwrite(helper,sizeof(int),2*(nPDFx*nPDFy*nPDFz),write_ptr);
    fclose(write_ptr);

    //write id file
    sprintf(idname, "../data/tracer-%.4E/sortedid.%05d", *ts, *myid);
    write_ptr = fopen(idname, "wb");
    fwrite(ssn_out,sizeof(int64_t),(*fill),write_ptr);
    fclose(write_ptr);

    //write global position file
    sprintf(xlocname, "../data/tracer-%.4E/sortedposition.%05d", *ts, *myid);
    write_ptr = fopen(xlocname, "wb");
    fwrite(xloc_out,sizeof(double),3*(*fill),write_ptr);
    fclose(write_ptr);

    //write tracer variable file
    sprintf(tracername, "../data/tracer-%.4E/sortedtracer.%05d", *ts, *myid);
    write_ptr = fopen(tracername, "wb");
    fwrite(sfield_out,sizeof(double),(*nfield)*(*fill),write_ptr);
    fclose(write_ptr);
*/



    // stop profiler
    profiler_pause_("new output particles");

    //free memory
    free(PDFindex);
    free(helper);
    free(sorter);
    free(ssn_out);
    free(xloc_out);
    free(sfield_out);
    return;
}
