/* This code is part of the tng binary trajectory format.
 *
 *  The high-level API of the TNG API is used where appropriate.
 *
 * Written by Magnus Lundborg
 * Copyright (c) 2012-2013, The GROMACS development team.
 * Check out http://www.gromacs.org for more information.
 *
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the Revised BSD License.
 *
 * This test program is written by Hugo MacDermott-Opeskin 2020
 */

// is called extern C if __cplusplus is defined, so no need to call it explicitly here //
#include "tng/tng_io.h"

#ifdef USE_STD_INTTYPES_H
#    include <inttypes.h>
#endif

#include <string>
#include <iostream>

int main(int argc, char* argv[])
{

    // parse command line args
    if (argc <= 1)
    {
        printf(" no file specified\n usage tng_io_high_level <tng file> [-v] [start=0] [end=5000] "
               "\n");
        exit(0);
    }
    if (argc > 5)
    {
        printf("too many arguments");
        exit(0);
    }

    // C function signature requires char*
    char* trajname = argv[1];

    // parse flag for verbose
    bool verbose = false;

    char* flag = argv[2];
    if (flag)
    {
        printf("Flag %s set \n", flag);
        if (strcmp(flag, "-v") != 0)
        {
            printf("\ninvalid additional flags\nusage tng_io_high_level <tng file> [-v] [start=0] "
                   "[end=5000] \n");
            exit(0);
        }
        else
        {
            printf("setting verbose\n");
            verbose = true;
        }
    }

    // set the frame range, default is start = 0 and end = 5000
    int64_t start = 0;
    int64_t end   = 5000;

    // parse frame range if given
    if (argc > 3)
    {
        start = std::strtol(argv[3], nullptr, 10);
        if (argc > 4)
        {
            end = std::strtol(argv[4], nullptr, 10);
        }
    }

    int64_t nfr = end - start + 1;

    printf("TNG IO READ TESTS\n");
    printf("-----------------\n\n");
    printf("Trajectory %s \n", trajname);
    printf("reading between %ld and %ld TNG frames (steps) with number of frames %ld \n", start, end, nfr);
    printf("-----------------\n\n");

    // values into which we put the trajectory data
    // C function signatures requires *float
    float * positions = nullptr, *box_shape = nullptr, *forces = nullptr, *velocities = nullptr;
    int64_t n_particles, n_frames, tot_n_steps, i, j, stride_length, n_strides;

    // read
    tng_trajectory_t traj;
    tng_util_trajectory_open(trajname, 'r', &traj);

    // num_particles
    if (tng_num_particles_get(traj, &n_particles) != TNG_SUCCESS)
    {
        printf("Test read particles: FAIL\n");
        tng_util_trajectory_close(&traj);
    }
    else
    {
        printf("Test read particle:  PASS\n");
        if (verbose)
        {
            printf("N particles  %ld \n", n_particles);
        }
    }
    printf("-------------------------\n\n");


    // frames
    if (tng_num_frames_get(traj, &tot_n_steps) != TNG_SUCCESS)
    {
        printf("Test read frames:    FAIL\n");
        tng_util_trajectory_close(&traj);
    }
    else
    {
        printf("Test read frames:    PASS\n");
        if (verbose)
        {
            printf("N frames %ld \n", tot_n_steps);
        }
    }
    printf("-------------------------\n\n");


    // box shape
    if (tng_util_box_shape_read_range(traj, start, end, &box_shape, &stride_length) != TNG_SUCCESS)
    {
        printf("Test read box:       FAIL\n");
    }
    else
    {
        printf("Test read box:       PASS\n");
        if (verbose)
        {
            // figure out the number of frames from the number of steps
            n_frames = (nfr - 1) / stride_length + 1;
            printf("number of frames with box data = %ld\n", n_frames);
            printf("box data stride length = %ld\n", stride_length);

            printf("Simulation box shape:\n");
            for (i = 0; i < n_frames; i++)
            {
                printf("idx %ld frame %ld ", i, start + i * stride_length);
                for (j = 0; j < 9; j++)
                {
                    printf("%f ", box_shape[9 * i + j]);
                }
                printf("\n");
            }
        }
    }
    printf("-------------------------\n\n");


    // positions
    if (tng_util_pos_read_range(traj, start, end, &positions, &stride_length) != TNG_SUCCESS)
    {
        printf("Test read positions: FAIL\n");
    }
    else
    {
        printf("Test read positions: PASS\n");
        if (verbose)
        {
            // figure out the number of frames from the number of steps
            n_frames = (nfr - 1) / stride_length + 1;

            printf("number of frames with position data = %ld  in frame range %ld to %ld \n",
                   n_frames, start, end);
            printf("position data stride length = %ld\n", stride_length);

            // print the coords
            for (i = 0; i < n_frames; i++)
            {
                printf("idx %ld frame %ld \n", i, start + i * stride_length);
                for (j = 0; j < n_particles; j++)
                {
                    printf("%ld %ld %f  ", i, j, positions[i * n_particles * 3 + j * 3]);
                    printf("%ld %ld %f  ", i, j, positions[i * n_particles * 3 + j * 3 + 1]);
                    printf("%ld %ld %f  %ld \n", i, j, positions[i * n_particles * 3 + j * 3 + 2],
                           i * n_particles * 3 + j * 3 + 2);
                }
            }
        }
    }
    printf("-------------------------\n\n");


    // forces
    if (tng_util_force_read(traj, &forces, &stride_length) != TNG_SUCCESS)
    {
        printf("Test read forces:    FAIL\n");
    }
    else
    {
        printf("Test read forces:    PASS\n");
    }
    printf("-------------------------\n\n");

    // velocities
    if (tng_util_vel_read(traj, &velocities, &stride_length) != TNG_SUCCESS)
    {
        printf("Test read vels:      FAIL\n");
    }
    else
    {
        printf("Test read vels:      PASS\n");
    }
    printf("-------------------------\n\n");


    // clean up the C memory

    if (positions)
    {
        printf("free positions\n");
        free(positions);
    }
    if (box_shape)
    {
        printf("free box\n");
        free(box_shape);
    }

    if (forces)
    {
        printf("free forces\n");
        free(forces);
    }

    if (velocities)
    {
        printf("free vels\n");
        free(velocities);
    }
    tng_util_trajectory_close(&traj);
    return (0);
}
