/*
 * This code is part of the tng binary trajectory format.
 *
 * Copyright (c) 2012-2013, The GROMACS development team.
 * Copyright (c) 2020, by the GROMACS development team.
 * TNG was orginally written by Magnus Lundborg, Daniel Spångberg and
 * Rossen Apostolov. The API is implemented mainly by Magnus Lundborg,
 * Daniel Spångberg and Anders Gärdenäs.
 *
 * Please see the AUTHORS file for more information.
 *
 * The TNG library is free software; you can redistribute it and/or
 * modify it under the terms of the Revised BSD License.
 *
 * To help us fund future development, we humbly ask that you cite
 * the research papers on the package.
 *
 * Check out http://www.gromacs.org for more information.
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


    printf("TNG IO READ TESTS\n");
    printf("-----------------\n\n");
    printf("Trajectory %s \n", trajname);
    printf("reading between %ld and %ld TNG frames (steps) specified on command line \n", start, end);
    printf("-----------------\n\n");

    // values into which we put the trajectory data
    // C function signatures requires *float
    float * positions = nullptr, *box_shape = nullptr, *forces = nullptr, *velocities = nullptr;
    int64_t n_particles, n_frames, tot_n_steps, i, j, stride_length, n_strides;
    tng_function_status stat;

    // read
    tng_trajectory_t traj;
    stat = tng_util_trajectory_open(trajname, 'r', &traj);

    // failure or critical is exit(1)
    if (stat != TNG_SUCCESS)
    {
        printf("Test Trajectory Open: FAIL\n");
        if (stat != TNG_FAILURE)
        {
            printf("Failure is TNG_CRITICAL\n");
        }
        tng_util_trajectory_close(&traj);
        exit(1);
    }
    else
    {
        printf("Test traj open:      PASS\n");
    }
    printf("-------------------------\n\n");

    // failure or critical is exit(1)
    stat = tng_num_particles_get(traj, &n_particles);
    if (stat != TNG_SUCCESS)
    {
        printf("Test read particles: FAIL\n");
        if (stat != TNG_FAILURE)
        {
            printf("Failure is TNG_CRITICAL\n");
        }
        tng_util_trajectory_close(&traj);
        exit(1);
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


    // failure or critical is exit(1)
    stat = tng_num_frames_get(traj, &tot_n_steps);
    if (stat != TNG_SUCCESS)
    {
        printf("Test read frames:    FAIL\n");
        if (stat != TNG_FAILURE)
        {
            printf("Failure is TNG_CRITICAL\n");
        }
        tng_util_trajectory_close(&traj);
        exit(1);
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

    // correct for if the number of frames is less than the default
    if (end > tot_n_steps - 1)
    {
        end = tot_n_steps - 1;
        if (verbose)
        {
            printf("Frame range clipped to shorter trajectory \nFrame range to actually read "
                   "between  start = %ld end = %ld \n\n",
                   start, end);
        }
    }

    int64_t nfr = end - start + 1;


    // critical is exit(1)
    stat = tng_util_box_shape_read_range(traj, start, end, &box_shape, &stride_length);
    if (stat != TNG_SUCCESS)
    {
        printf("Test read box:       FAIL\n");
        if (stat != TNG_FAILURE)
        {
            printf("Failure is TNG_CRITICAL\n");
            tng_util_trajectory_close(&traj);
            exit(1);
        }
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


    // critical is exit(1)
    stat = tng_util_pos_read_range(traj, start, end, &positions, &stride_length);
    if (stat != TNG_SUCCESS)
    {
        printf("Test read positions: FAIL\n");
        if (stat != TNG_FAILURE)
        {
            printf("Failure is TNG_CRITICAL\n");
            tng_util_trajectory_close(&traj);
            exit(1);
        }
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


    // critical is exit(1)
    stat = tng_util_force_read_range(traj, start, end, &forces, &stride_length);
    if (stat != TNG_SUCCESS)
    {
        printf("Test read forces:    FAIL\n");
        if (stat != TNG_FAILURE)
        {
            printf("Failure is TNG_CRITICAL\n");
            tng_util_trajectory_close(&traj);
            exit(1);
        }
    }
    else
    {
        printf("Test read forces:    PASS\n");
    }
    printf("-------------------------\n\n");

    // velocities
    stat = tng_util_vel_read_range(traj, start, end, &velocities, &stride_length);
    if (stat != TNG_SUCCESS)
    {
        printf("Test read vels:      FAIL\n");
        if (stat != TNG_FAILURE)
        {
            printf("Failure is TNG_CRITICAL\n");
            tng_util_trajectory_close(&traj);
            exit(1);
        }
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
