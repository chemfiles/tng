/* This code is part of the tng binary trajectory format.
 *
 *  The high-level API of the TNG API is used where appropriate.
 *
 *                      VERSION 1.0
 *
 * Written by Magnus Lundborg
 * Copyright (c) 2012, The GROMACS development team.
 * check out http://www.gromacs.org for more information.
 *
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 */

#include <stdlib.h>
#include <stdio.h>
#include <tng_io.h>

int main(int argc, char **argv)
{
    tng_trajectory_t traj;
    /* Assume that the data is stored as floats. The data is placed in 1-D
     * arrays */
    float *positions = 0, *box_shape = 0;
    int64_t n_particles, n_frames, tot_n_frames, stride_length, i, j;
    // Set a default frame range
    int64_t first_frame = 0, last_frame = 5000, n_strides;

    if(argc <= 1)
    {
        printf("No file specified\n");
        printf("Usage:\n");
        printf("tng_io_read_pos <tng_file> "
               "[first_frame = %"PRId64"] [last_frame = %"PRId64"]\n",
               first_frame, last_frame);
        exit(1);
    }

    // A reference must be passed to allocate memory
    tng_util_trajectory_open(argv[1], 'r', &traj);

    if(argc >= 3)
    {
        first_frame = strtoll(argv[2], 0, 10);
        if(argc >= 4)
        {
            last_frame = strtoll(argv[3], 0, 10);
        }
    }

    if(tng_num_frames_get(traj, &tot_n_frames) != TNG_SUCCESS)
    {
        printf("Cannot determine the number of frames in the file\n");
        tng_util_trajectory_close(&traj);
        exit(1);
    }

    if(tng_num_particles_get(traj, &n_particles) != TNG_SUCCESS)
    {
        printf("Cannot determine the number of particles in the file\n");
        tng_util_trajectory_close(&traj);
        exit(1);
    }

    printf("%"PRId64" frames in file\n", tot_n_frames);

    if(last_frame > tot_n_frames - 1)
    {
        last_frame = tot_n_frames - 1;
    }


    n_frames = last_frame - first_frame + 1;

    if(tng_util_box_shape_read(traj, &box_shape, &stride_length) ==
       TNG_SUCCESS)
    {
        printf("Simulation box shape: ");
        for(i=0; i < 9; i++)
        {
            printf("%f ", box_shape[i]);
        }
        printf("\n");
    }
    else
    {
        printf("Simulation box shape not set in the file (or could not be read)\n");
    }

    n_strides = (n_frames % stride_length) ? n_frames / stride_length + 1:
                n_frames / stride_length;

    // Get the positions of all particles in the requested frame range.
    // The positions are stored in the positions array.
    // N.B. No proper error checks.
    if(tng_util_pos_read_range(traj, 0, last_frame, &positions, &stride_length)
       == TNG_SUCCESS)
    {
        // Print the positions of the wanted particle (zero based)
        for(i=0; i < n_strides; i++)
        {
            printf("\nFrame %"PRId64":\n", first_frame + i*stride_length);
            for(j=0; j < n_particles; j++)
            {
                printf("Atom nr: %"PRId64"", j);
                printf("\t%f", positions[i*n_particles*3+j*3]);
                printf("\t%f", positions[i*n_particles*3+j*3+1]);
                printf("\t%f", positions[i*n_particles*3+j*3+2]);
                printf("\n");
            }
        }
    }
    else
    {
        printf("Cannot read positions\n");
    }

    // Free memory
    if(positions)
    {
        free(positions);
    }
    if(box_shape)
    {
        free(box_shape);
    }
    tng_util_trajectory_close(&traj);

    return(0);
}
