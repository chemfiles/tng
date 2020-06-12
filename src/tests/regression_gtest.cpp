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


#include <iostream>
#include <fstream>

// is called extern C if __cplusplus is defined, so no need to call it explicitly here //
#include "tng/tng_io.h"
#include "gtest/gtest.h"

#ifdef USE_STD_INTTYPES_H
#    include <inttypes.h>
#endif

// space for helper functions


// linked with gtest_main which returns 0 if all tests pass

// build class for each trajectory file
class ArgonCompressedTest : public ::testing::Test
{
protected:
    tng_trajectory_t    traj;
    tng_function_status stat_init;
    const char*         filename = "./example_files/argon_npt_compressed.tng";

    ArgonCompressedTest() { stat_init = tng_util_trajectory_open(filename, 'r', &traj); }

    ~ArgonCompressedTest() override { tng_util_trajectory_close(&traj); }
};


TEST_F(ArgonCompressedTest, TrajectoryOpen)
{

    EXPECT_EQ(stat_init, TNG_SUCCESS);
}

TEST_F(ArgonCompressedTest, NumParticles)
{
    int64_t             n_particles;
    tng_function_status read_stat;
    read_stat = tng_num_particles_get(traj, &n_particles);
    EXPECT_EQ(read_stat, TNG_SUCCESS);
    EXPECT_EQ(n_particles, 1000);
}

TEST_F(ArgonCompressedTest, NumFrames)
{
    int64_t             n_frames;
    tng_function_status read_stat;
    read_stat = tng_num_frames_get(traj, &n_frames);
    EXPECT_EQ(read_stat, TNG_SUCCESS);
    EXPECT_EQ(n_frames, 500001);
}

TEST_F(ArgonCompressedTest, BoxShapeRead)
{
    float*              box_shape = nullptr;
    int64_t             stride_length;
    tng_function_status read_stat;
    read_stat = tng_util_box_shape_read(traj, &box_shape, &stride_length);
    EXPECT_EQ(read_stat, TNG_SUCCESS);
    EXPECT_EQ(stride_length, 5000);
}


// incomplete
TEST_F(ArgonCompressedTest, BoxShapeValues)
{
    float*              box_shape = nullptr;
    int64_t             stride_length, i;
    tng_function_status read_stat;
    read_stat = tng_util_box_shape_read(traj, &box_shape, &stride_length);
    EXPECT_EQ(read_stat, TNG_SUCCESS);
    // box_shape frame 0
    // clang-format off
    const std::vector<float> frame_0  = {
        3.60140, 0.00000, 0.000000, 0.000000, 3.60140, 0.000000, 0.000000, 0.000000, 3.60140,
    };
    // clang-format on
    ASSERT_EQ(frame_0.size(), 9);
    // compare element by element
    for (i = 0; i < 9; i++)
    {
        EXPECT_FLOAT_EQ(box_shape[i], frame_0[i]);
    }
    const std::vector<float> frame_100 = {
        3.589650, 0.000000, 0.000000, 0.000000, 3.589650, 0.000000, 0.000000, 0.000000, 3.589650,
    };
    ASSERT_EQ(frame_100.size(), 9);
    for (i = 0; i < 9; i++)
    {
        EXPECT_FLOAT_EQ(box_shape[909 - 9 + i], frame_100[i]);
    }
}

TEST_F(ArgonCompressedTest, PositionPartialRead)
{
    float*              positions = nullptr;
    int64_t             stride_length;
    tng_function_status read_stat;
    read_stat = tng_util_pos_read_range(traj, 0, 5001, &positions, &stride_length);
    EXPECT_EQ(read_stat, TNG_SUCCESS);
    EXPECT_EQ(stride_length, 5000);
}

TEST_F(ArgonCompressedTest, PositionPartialReadIrregular)
{
    float*              positions = nullptr;
    int64_t             stride_length;
    tng_function_status read_stat;
    read_stat = tng_util_pos_read_range(traj, 7777, 18888, &positions, &stride_length);
    EXPECT_EQ(read_stat, TNG_SUCCESS);
    EXPECT_EQ(stride_length, 5000);
}

TEST_F(ArgonCompressedTest, PositionRead)
{
    float*              positions = nullptr;
    int64_t             stride_length;
    tng_function_status read_stat;
    read_stat = tng_util_pos_read(traj, &positions, &stride_length);
    EXPECT_EQ(read_stat, TNG_SUCCESS);
    EXPECT_EQ(stride_length, 5000);
}
// incomplete
TEST_F(ArgonCompressedTest, PositionValues)
{
    float*              positions = nullptr;
    int64_t             start, end, stride_length, i;
    tng_function_status read_stat;

    // xyz first 10 atoms frame 0
    // clang-format off
    const std::vector<float> frame_0_first_10_values  = {
    2.53300e+00,  1.24400e+00,  3.50600e+00,
    8.30000e-01,  2.54400e+00,  3.44800e+00,
    1.09100e+00,  1.10000e-01,  3.12900e+00,
    2.45500e+00,  5.00000e-03,  3.01200e+00,
    2.71400e+00,  1.35300e+00,  5.53000e-01,
    3.05100e+00,  2.89300e+00,  2.69100e+00,
    1.42200e+00,  2.77000e+00,  1.46000e-01,
    2.22300e+00,  1.21100e+00,  3.26800e+00,
    2.81100e+00,  2.78900e+00,  2.38500e+00,
    4.87000e-01,  1.15900e+00,  1.17100e+00
    };
    // clang-format on
    ASSERT_EQ(frame_0_first_10_values.size(), 30);
    read_stat = tng_util_pos_read(traj, &positions, &stride_length);
    EXPECT_EQ(read_stat, TNG_SUCCESS);
    // compare element by element
    for (i = 0; i < 30; i++)
    {
        EXPECT_FLOAT_EQ(positions[i], frame_0_first_10_values[i]);
    }


    // xyz last 10 atoms frame 100
    // clang-format off
    const std::vector<float> frame_100_last_10_values  = {
    0.776000,  1.196000,   0.773000,  
    0.627000,  0.334000,   2.049000,  
    0.609000,  3.463000,   0.257000,  
    3.020000,  3.184000,   2.976000,   
    2.647000,  0.774000,   1.815000,   
    0.156000,  1.283000,   3.281000,  
    0.658000,  3.033000,   2.908000,  
    2.085000,  3.551000,   1.436000,  
    0.156000,  3.502000,   0.314000,  
    1.289000,  0.998000,   1.645000  
    };
    // clang-format on
    EXPECT_EQ(frame_100_last_10_values.size(), 30);
    // compare element by element
    for (i = 0; i < 30; i++)
    {
        EXPECT_FLOAT_EQ(positions[303000 - 30 + i], frame_100_last_10_values[i]);
    }
}

TEST_F(ArgonCompressedTest, ForceRead)
{
    float*              forces = nullptr;
    int64_t             stride_length;
    tng_function_status read_stat;
    read_stat = tng_util_force_read(traj, &forces, &stride_length);
    // no forces expexted in this file
    EXPECT_EQ(read_stat, TNG_FAILURE);
}

TEST_F(ArgonCompressedTest, VelRead)
{
    float*              velocities = nullptr;
    int64_t             stride_length;
    tng_function_status read_stat;
    read_stat = tng_util_vel_read(traj, &velocities, &stride_length);
    // no velocities expected in this file
    EXPECT_EQ(read_stat, TNG_FAILURE);
}

TEST_F(ArgonCompressedTest, NumMoleculeTypes)
{
    int64_t             count;
    tng_function_status read_stat;
    read_stat = tng_num_molecule_types_get(traj, &count);
    EXPECT_EQ(read_stat, TNG_SUCCESS);
    EXPECT_EQ(count, 1);
}

TEST_F(ArgonCompressedTest, NumMolecules)
{
    int64_t             count;
    tng_function_status read_stat;
    read_stat = tng_num_molecules_get(traj, &count);
    EXPECT_EQ(read_stat, TNG_SUCCESS);
    EXPECT_EQ(count, 1000);
}

TEST_F(ArgonCompressedTest, MoleculeFind)
{
    int64_t             count;
    tng_molecule_t      molecule;
    tng_function_status read_stat;
    read_stat = tng_molecule_find(traj, "Argon", -1, &molecule);
    EXPECT_EQ(read_stat, TNG_SUCCESS);
}


class TNGExampleTest : public ::testing::Test
{
protected:
    tng_trajectory_t    traj;
    tng_function_status stat_init;
    const char*         filename = "./example_files/tng_example.tng";

    TNGExampleTest() { stat_init = tng_util_trajectory_open(filename, 'r', &traj); }

    ~TNGExampleTest() override { tng_util_trajectory_close(&traj); }
};


TEST_F(TNGExampleTest, TrajectoryOpen)
{

    EXPECT_EQ(stat_init, TNG_SUCCESS);
}

TEST_F(TNGExampleTest, NumParticles)
{
    int64_t             n_particles;
    tng_function_status read_stat;
    read_stat = tng_num_particles_get(traj, &n_particles);
    EXPECT_EQ(read_stat, TNG_SUCCESS);
    EXPECT_EQ(n_particles, 15);
}

TEST_F(TNGExampleTest, NumFrames)
{
    int64_t             n_frames;
    tng_function_status read_stat;
    read_stat = tng_num_frames_get(traj, &n_frames);
    EXPECT_EQ(read_stat, TNG_SUCCESS);
    EXPECT_EQ(n_frames, 10);
}

TEST_F(TNGExampleTest, BoxShapeRead)
{
    float*              box_shape = nullptr;
    int64_t             stride_length;
    tng_function_status read_stat;
    read_stat = tng_util_box_shape_read(traj, &box_shape, &stride_length);
    EXPECT_EQ(read_stat, TNG_FAILURE);
}

TEST_F(TNGExampleTest, BoxReadThenPositionPartialRead)
{
    float*              positions = nullptr;
    float*              box_shape = nullptr;
    int64_t             stride_length;
    tng_function_status read_stat;
    EXPECT_EQ(stat_init, TNG_SUCCESS);
    read_stat = tng_util_box_shape_read(traj, &box_shape, &stride_length);
    // normal box not expected
    EXPECT_EQ(read_stat, TNG_FAILURE);
    read_stat = tng_util_pos_read_range(traj, 0, 2, &positions, &stride_length);
    EXPECT_EQ(read_stat, TNG_SUCCESS);
    EXPECT_EQ(stride_length, 1);
}


TEST_F(TNGExampleTest, PositionPartialRead)
{
    float*              positions = nullptr;
    int64_t             stride_length;
    tng_function_status read_stat;
    EXPECT_EQ(stat_init, TNG_SUCCESS);
    read_stat = tng_util_pos_read_range(traj, 0, 2, &positions, &stride_length);
    EXPECT_EQ(read_stat, TNG_SUCCESS);
    EXPECT_EQ(stride_length, 1);
}


TEST_F(TNGExampleTest, PositionRead)
{
    float*              positions = nullptr;
    int64_t             stride_length;
    tng_function_status read_stat;
    read_stat = tng_util_pos_read(traj, &positions, &stride_length);
    EXPECT_EQ(read_stat, TNG_SUCCESS);
    EXPECT_EQ(stride_length, 1);
}
// incomplete
TEST_F(TNGExampleTest, PositionValues)
{
    float*              positions = nullptr;
    int64_t             stride_length, i;
    tng_function_status read_stat;
    read_stat = tng_util_pos_read(traj, &positions, &stride_length);
    EXPECT_EQ(read_stat, TNG_SUCCESS);

    // xyz first 10 atoms frame 0
    // clang-format off
    const std::vector<float> frame_0_first_10_values  = {
    1.000000,   1.000000,   1.000000,   
    2.000000,   2.000000,   2.000000,  
    3.000000,   3.000000,   3.000000,  
    11.000000,  11.000000,  11.000000,  
    12.000000,  12.000000,  12.000000,  
    13.000000,  13.000000,  13.000000,  
    21.000000,  21.000000,  21.000000,  
    22.000000,  22.000000,  22.000000,  
    23.000000,  23.000000,  23.000000,  
    8.250000,   33.000000,  33.000000,  
    };
    // clang-format on
    ASSERT_EQ(frame_0_first_10_values.size(), 30);
    // compare element by element
    for (i = 0; i < 30; i++)
    {
        EXPECT_FLOAT_EQ(positions[i], frame_0_first_10_values[i]);
    }
    // clang-format off
    const std::vector<float> frame_9_last_10_values  = {
    13.000000,  13.000000,  13.000000,  
    21.000000,  21.000000,  21.000000,
    22.000000,  22.000000,  22.000000,  
    23.000000,  23.000000,  23.000000,  
    8.250000,   33.000000,  33.000000,  
    8.250000,   34.000000,  33.000000, 
    8.500000,   33.000000,  34.000000, 
    50.000000,  50.000000,  50.000000, 
    51.000000,  51.000000,  51.000000, 
    100.000000, 100.000000, 100.000000,
    };
    // clang-format on
    ASSERT_EQ(frame_0_first_10_values.size(), 30);
    // compare element by element
    for (i = 0; i < 30; i++)
    {
        EXPECT_FLOAT_EQ(positions[450 - 30 + i], frame_9_last_10_values[i]);
    }
}

TEST_F(TNGExampleTest, ForceRead)
{
    float*              forces = nullptr;
    int64_t             stride_length;
    tng_function_status read_stat;
    read_stat = tng_util_force_read(traj, &forces, &stride_length);
    // no forces expected in this file
    EXPECT_EQ(read_stat, TNG_FAILURE);
}

TEST_F(TNGExampleTest, VelRead)
{
    float*              velocities = nullptr;
    int64_t             stride_length;
    tng_function_status read_stat;
    read_stat = tng_util_vel_read(traj, &velocities, &stride_length);
    // no velocities expected in this file
    EXPECT_EQ(read_stat, TNG_FAILURE);
}

TEST_F(TNGExampleTest, NumMoleculeTypes)
{
    int64_t             count;
    tng_function_status read_stat;
    read_stat = tng_num_molecule_types_get(traj, &count);
    EXPECT_EQ(read_stat, TNG_SUCCESS);
    EXPECT_EQ(count, 1);
}


TEST_F(TNGExampleTest, NumMolecules)
{
    int64_t             count;
    tng_function_status read_stat;
    read_stat = tng_num_molecules_get(traj, &count);
    EXPECT_EQ(read_stat, TNG_SUCCESS);
    // 5 water molecules with 3 particles each
    EXPECT_EQ(count, 5);
}

TEST_F(TNGExampleTest, MoleculeFind)
{
    int64_t             count;
    tng_molecule_t      molecule;
    tng_function_status read_stat;
    read_stat = tng_molecule_find(traj, "water", -1, &molecule);
    EXPECT_EQ(read_stat, TNG_SUCCESS);
}
