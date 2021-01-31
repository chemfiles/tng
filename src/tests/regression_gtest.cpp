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

void free_float_data_if_present(float* data)
{
    if (data)
    {
        free(data);
    }
}

// linked with gtest_main which returns 0 if all tests pass

// build class for each trajectory file
class ArgonCompressedTest : public ::testing::Test
{
protected:
    tng_trajectory_t    traj;
    tng_function_status stat_init;
    float*              positions  = nullptr;
    float*              box_shape  = nullptr;
    float*              forces     = nullptr;
    float*              velocities = nullptr;
    const char*         filename   = "./example_files/argon_npt_compressed.tng";

    ArgonCompressedTest() { stat_init = tng_util_trajectory_open(filename, 'r', &traj); }

    ~ArgonCompressedTest() override
    {
        tng_util_trajectory_close(&traj);
        free_float_data_if_present(positions);
        free_float_data_if_present(box_shape);
        free_float_data_if_present(forces);
        free_float_data_if_present(velocities);
    }
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
    int64_t             stride_length;
    tng_function_status read_stat;
    read_stat = tng_util_box_shape_read(traj, &box_shape, &stride_length);
    EXPECT_EQ(read_stat, TNG_SUCCESS);
    EXPECT_EQ(stride_length, 5000);
}

TEST_F(ArgonCompressedTest, BoxShapeValues)
{
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
    int64_t             stride_length;
    tng_function_status read_stat;
    read_stat = tng_util_pos_read_range(traj, 0, 5000, &positions, &stride_length);
    EXPECT_EQ(read_stat, TNG_SUCCESS);
    EXPECT_EQ(stride_length, 5000);
}

TEST_F(ArgonCompressedTest, PositionPartialReadInvalidRange)
{
    int64_t             stride_length;
    tng_function_status read_stat;
    read_stat = tng_util_pos_read_range(traj, 1, 1, &positions, &stride_length);
    EXPECT_EQ(read_stat, TNG_FAILURE);
    EXPECT_EQ(stride_length, 5000);
}

TEST_F(ArgonCompressedTest, PositionPartialValuesFrm0)
{
    int64_t             stride_length, i;
    tng_function_status read_stat;
    read_stat = tng_util_pos_read_range(traj, 0, 0, &positions, &stride_length);
    EXPECT_EQ(read_stat, TNG_SUCCESS);
    EXPECT_EQ(stride_length, 5000);
    // xyz first 10 atoms frame 0
    // gmx dump frame 0
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
    4.87000e-01,  1.15900e+00,  1.17100e+00,
    };
    // clang-format on
    ASSERT_EQ(frame_0_first_10_values.size(), 30);
    // compare element by element
    for (i = 0; i < 30; i++)
    {
        EXPECT_FLOAT_EQ(positions[i], frame_0_first_10_values[i]);
    }
}


TEST_F(ArgonCompressedTest, PositionPartialValuesFrm1)
{
    int64_t             stride_length, i;
    tng_function_status read_stat;
    read_stat = tng_util_pos_read_range(traj, 5000, 5000, &positions, &stride_length);
    EXPECT_EQ(read_stat, TNG_SUCCESS);
    EXPECT_EQ(stride_length, 5000);
    // xyz first 10 atoms frame 1
    // gmx dump frame 1
    // clang-format off
    const std::vector<float> frame_1_first_10_values  = {
    2.52400e+00,  1.18600e+00,  2.33000e-01,
    9.01000e-01,  2.77300e+00,  3.13500e+00,
    1.68400e+00,  2.14000e-01,  3.17900e+00,
    2.19300e+00,  3.37400e+00,  2.90700e+00,
    2.89400e+00,  1.50600e+00,  4.46000e-01,
    2.79600e+00,  2.87600e+00,  2.54400e+00,
    1.37600e+00,  3.06700e+00,  1.78000e-01,
    2.17900e+00,  9.04000e-01,  3.04800e+00,
    2.68400e+00,  2.53800e+00,  2.47300e+00,
    3.78000e-01,  1.41000e+00,  8.46000e-01,
    };
    // clang-format on
    ASSERT_EQ(frame_1_first_10_values.size(), 30);
    // compare element by element
    for (i = 0; i < 30; i++)
    {
        EXPECT_FLOAT_EQ(positions[i], frame_1_first_10_values[i]);
    }
}

TEST_F(ArgonCompressedTest, PositionPartialReadIrregular)
{
    int64_t             stride_length;
    tng_function_status read_stat;
    read_stat = tng_util_pos_read_range(traj, 7777, 18888, &positions, &stride_length);
    EXPECT_EQ(read_stat, TNG_SUCCESS);
    EXPECT_EQ(stride_length, 5000);
}

TEST_F(ArgonCompressedTest, PositionRead)
{
    int64_t             stride_length;
    tng_function_status read_stat;
    read_stat = tng_util_pos_read(traj, &positions, &stride_length);
    EXPECT_EQ(read_stat, TNG_SUCCESS);
    EXPECT_EQ(stride_length, 5000);
}

TEST_F(ArgonCompressedTest, PositionValues)
{
    int64_t             start, end, stride_length, i;
    tng_function_status read_stat;
    read_stat = tng_util_pos_read(traj, &positions, &stride_length);
    EXPECT_EQ(read_stat, TNG_SUCCESS);
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
    4.87000e-01,  1.15900e+00,  1.17100e+00,
    };
    // clang-format on
    ASSERT_EQ(frame_0_first_10_values.size(), 30);
    // compare element by element
    for (i = 0; i < 30; i++)
    {
        EXPECT_FLOAT_EQ(positions[i], frame_0_first_10_values[i]);
    }


    // xyz last 10 atoms frame 100
    // clang-format off
    const std::vector<float> frame_100_last_10_values  = {
    7.76000e-01,  1.19600e+00,  7.73000e-01, 
    6.27000e-01,  3.34000e-01,  2.04900e+00, 
    6.09000e-01,  3.46300e+00,  2.57000e-01, 
    3.02000e+00,  3.18400e+00,  2.97600e+00,  
    2.64700e+00,  7.74000e-01,  1.81500e+00,  
    1.56000e-01,  1.28300e+00,  3.28100e+00, 
    6.58000e-01,  3.03300e+00,  2.90800e+00, 
    2.08500e+00,  3.55100e+00,  1.43600e+00, 
    1.56000e-01,  3.50200e+00,  3.14000e-01, 
    1.28900e+00,  9.98000e-01,  1.64500e+00,
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
    int64_t             stride_length;
    tng_function_status read_stat;
    read_stat = tng_util_force_read(traj, &forces, &stride_length);
    // no forces expexted in this file
    EXPECT_EQ(read_stat, TNG_FAILURE);
}

TEST_F(ArgonCompressedTest, VelRead)
{
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
    float*              positions  = nullptr;
    float*              box_shape  = nullptr;
    float*              forces     = nullptr;
    float*              velocities = nullptr;
    const char*         filename   = "./example_files/tng_example.tng";

    TNGExampleTest() { stat_init = tng_util_trajectory_open(filename, 'r', &traj); }

    ~TNGExampleTest() override
    {
        tng_util_trajectory_close(&traj);
        free_float_data_if_present(positions);
        free_float_data_if_present(box_shape);
        free_float_data_if_present(forces);
        free_float_data_if_present(velocities);
    }
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
    int64_t             stride_length;
    tng_function_status read_stat;
    read_stat = tng_util_box_shape_read(traj, &box_shape, &stride_length);
    EXPECT_EQ(read_stat, TNG_FAILURE);
}

TEST_F(TNGExampleTest, BoxReadThenPositionPartialRead)
{
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
    int64_t             stride_length;
    tng_function_status read_stat;
    EXPECT_EQ(stat_init, TNG_SUCCESS);
    read_stat = tng_util_pos_read_range(traj, 0, 2, &positions, &stride_length);
    EXPECT_EQ(read_stat, TNG_SUCCESS);
    EXPECT_EQ(stride_length, 1);
}


TEST_F(TNGExampleTest, PositionRead)
{
    int64_t             stride_length;
    tng_function_status read_stat;
    read_stat = tng_util_pos_read(traj, &positions, &stride_length);
    EXPECT_EQ(read_stat, TNG_SUCCESS);
    EXPECT_EQ(stride_length, 1);
}

TEST_F(TNGExampleTest, PositionValues)
{
    int64_t             stride_length, i;
    tng_function_status read_stat;
    read_stat = tng_util_pos_read(traj, &positions, &stride_length);
    EXPECT_EQ(read_stat, TNG_SUCCESS);

    // xyz first 10 atoms frame 0
    // gmx dump frame 0
    // clang-format off
    const std::vector<float> frame_0_first_10_values  = {
    1.00000e+00,  1.00000e+00,  1.00000e+00,
    2.00000e+00,  2.00000e+00,  2.00000e+00,
    3.00000e+00,  3.00000e+00,  3.00000e+00,
    1.10000e+01,  1.10000e+01,  1.10000e+01,
    1.20000e+01,  1.20000e+01,  1.20000e+01,
    1.30000e+01,  1.30000e+01,  1.30000e+01,
    2.10000e+01,  2.10000e+01,  2.10000e+01,
    2.20000e+01,  2.20000e+01,  2.20000e+01,
    2.30000e+01,  2.30000e+01,  2.30000e+01,
    8.25000e+00,  3.30000e+01,  3.30000e+01,
    };
    // clang-format on
    ASSERT_EQ(frame_0_first_10_values.size(), 30);
    // compare element by element
    for (i = 0; i < 30; i++)
    {
        EXPECT_FLOAT_EQ(positions[i], frame_0_first_10_values[i]);
    }
    // xyz first 10 atoms frame 9
    // gmx dump frame 9
    // clang-format off
    const std::vector<float> frame_9_last_10_values  = {
    1.30000e+01,  1.30000e+01,  1.30000e+01, 
    2.10000e+01,  2.10000e+01,  2.10000e+01,
    2.20000e+01,  2.20000e+01,  2.20000e+01, 
    2.30000e+01,  2.30000e+01,  2.30000e+01, 
    8.25000e+00,  3.30000e+01,  3.30000e+01, 
    8.25000e+00,  3.40000e+01,  3.30000e+01,
    8.50000e+00,  3.30000e+01,  3.40000e+01,
    5.00000e+01,  5.00000e+01,  5.00000e+01,
    5.10000e+01,  5.10000e+01,  5.10000e+01,
    1.00000e+02,  1.00000e+02,  1.00000e+02,
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
    int64_t             stride_length;
    tng_function_status read_stat;
    read_stat = tng_util_force_read(traj, &forces, &stride_length);
    // no forces expected in this file
    EXPECT_EQ(read_stat, TNG_FAILURE);
}

TEST_F(TNGExampleTest, VelRead)
{
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


class WaterTrjconvTest : public ::testing::Test
{
protected:
    tng_trajectory_t    traj;
    tng_function_status stat_init;
    float*              positions  = nullptr;
    float*              box_shape  = nullptr;
    float*              forces     = nullptr;
    float*              velocities = nullptr;
    const char*         filename   = "./example_files/water_npt_compressed_trjconv.tng";

    WaterTrjconvTest() { stat_init = tng_util_trajectory_open(filename, 'r', &traj); }

    ~WaterTrjconvTest() override
    {
        tng_util_trajectory_close(&traj);
        free_float_data_if_present(positions);
        free_float_data_if_present(box_shape);
        free_float_data_if_present(forces);
        free_float_data_if_present(velocities);
    }
};

TEST_F(WaterTrjconvTest, TrajectoryOpen)
{

    EXPECT_EQ(stat_init, TNG_SUCCESS);
}

TEST_F(WaterTrjconvTest, NumParticles)
{
    int64_t             n_particles;
    tng_function_status read_stat;
    read_stat = tng_num_particles_get(traj, &n_particles);
    EXPECT_EQ(read_stat, TNG_SUCCESS);
    EXPECT_EQ(n_particles, 2700);
}

TEST_F(WaterTrjconvTest, NumFrames)
{
    int64_t             n_frames;
    tng_function_status read_stat;
    read_stat = tng_num_frames_get(traj, &n_frames);
    EXPECT_EQ(read_stat, TNG_SUCCESS);
    EXPECT_EQ(n_frames, 500001);
}

TEST_F(WaterTrjconvTest, BoxShapeRead)
{
    int64_t             stride_length;
    tng_function_status read_stat;
    read_stat = tng_util_box_shape_read(traj, &box_shape, &stride_length);
    EXPECT_EQ(read_stat, TNG_SUCCESS);
    EXPECT_EQ(stride_length, 5000);
}


TEST_F(WaterTrjconvTest, BoxShapeValues)
{
    int64_t             stride_length, i;
    tng_function_status read_stat;
    read_stat = tng_util_box_shape_read(traj, &box_shape, &stride_length);
    EXPECT_EQ(read_stat, TNG_SUCCESS);
    const std::vector<float> frame_0 = {
        3.01125e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 3.01125e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 3.01125e+00,
    };
    ASSERT_EQ(frame_0.size(), 9);
    // compare element by element
    for (i = 0; i < 9; i++)
    {
        EXPECT_FLOAT_EQ(box_shape[i], frame_0[i]);
    }
    const std::vector<float> frame_100 = {
        2.870210, 0.000000, 0.000000, 0.000000, 2.870210, 0.000000, 0.000000, 0.000000, 2.870210,
    };
    ASSERT_EQ(frame_100.size(), 9);
    for (i = 0; i < 9; i++)
    {
        EXPECT_FLOAT_EQ(box_shape[909 - 9 + i], frame_100[i]);
    }
}

TEST_F(WaterTrjconvTest, PositionRead)
{
    int64_t             stride_length;
    tng_function_status read_stat;
    read_stat = tng_util_pos_read(traj, &positions, &stride_length);
    EXPECT_EQ(read_stat, TNG_SUCCESS);
    EXPECT_EQ(stride_length, 5000);
}

TEST_F(WaterTrjconvTest, PositionValues)
{
    int64_t             start, end, stride_length, i;
    tng_function_status read_stat;
    read_stat = tng_util_pos_read(traj, &positions, &stride_length);
    EXPECT_EQ(read_stat, TNG_SUCCESS);
    // xyz first 10 atoms frame 0
    // gmx dump frame 100
    // clang-format off
    const std::vector<float> frame_0_first_10_values  = {
    7.43000e-01,  2.42200e+00,  2.25100e+00,
    7.29000e-01,  2.48600e+00,  2.32000e+00,
    6.60000e-01,  2.37500e+00,  2.24500e+00,
    9.44000e-01,  1.43200e+00,  1.51800e+00,
    1.02100e+00,  1.48600e+00,  1.50000e+00,
    8.76000e-01,  1.46800e+00,  1.46000e+00,
    2.55700e+00,  2.11600e+00,  1.38800e+00,
    2.64500e+00,  2.14600e+00,  1.41200e+00,
    2.50000e+00,  2.15500e+00,  1.45400e+00,
    1.04500e+00,  2.51300e+00,  2.47000e-01,
    };
    // clang-format on
    ASSERT_EQ(frame_0_first_10_values.size(), 30);
    // compare element by element
    for (i = 0; i < 30; i++)
    {
        EXPECT_FLOAT_EQ(positions[i], frame_0_first_10_values[i]);
    }
    // xyz last 10 atoms frame 100
    // gmx dump frame 100
    // clang-format off
    const std::vector<float> frame_100_last_10_values  = {
    8.56000e-01,  2.24800e+00,  2.79100e+00,  
    9.53000e-01,  7.58000e-01,  8.59000e-01,  
    9.24000e-01,  7.24000e-01,  7.74000e-01,  
    8.73000e-01,  7.67000e-01,  9.10000e-01,  
    5.08000e-01,  2.31100e+00,  1.85000e-01,  
    5.25000e-01,  2.32700e+00,  9.30000e-02,  
    5.14000e-01,  2.21600e+00,  1.95000e-01,  
    6.67000e-01,  2.15500e+00,  1.65700e+00,  
    5.84000e-01,  2.14800e+00,  1.61000e+00,  
    7.18000e-01,  2.21700e+00,  1.60500e+00, 
    };
    // clang-format on
    ASSERT_EQ(frame_100_last_10_values.size(), 30);
    // compare element by element
    for (i = 0; i < 30; i++)
    {
        EXPECT_FLOAT_EQ(positions[818100 - 30 + i], frame_100_last_10_values[i]);
    }
}

TEST_F(WaterTrjconvTest, ForceRead)
{
    int64_t             stride_length;
    tng_function_status read_stat;
    read_stat = tng_util_force_read(traj, &forces, &stride_length);
    // no forces expexted in this file
    EXPECT_EQ(read_stat, TNG_FAILURE);
}

TEST_F(WaterTrjconvTest, VelRead)
{
    int64_t             stride_length;
    tng_function_status read_stat;
    read_stat = tng_util_vel_read(traj, &velocities, &stride_length);
    // no velocities expected in this file
    EXPECT_EQ(read_stat, TNG_FAILURE);
}

class WaterVelsForcesTest : public ::testing::Test
{
protected:
    tng_trajectory_t    traj;
    tng_function_status stat_init;
    float*              positions  = nullptr;
    float*              box_shape  = nullptr;
    float*              forces     = nullptr;
    float*              velocities = nullptr;
    const char*         filename   = "./example_files/water_uncompressed_vels_forces.tng";

    WaterVelsForcesTest() { stat_init = tng_util_trajectory_open(filename, 'r', &traj); }

    ~WaterVelsForcesTest() override
    {
        tng_util_trajectory_close(&traj);
        free_float_data_if_present(positions);
        free_float_data_if_present(box_shape);
        free_float_data_if_present(forces);
        free_float_data_if_present(velocities);
    }
};

TEST_F(WaterVelsForcesTest, TrajectoryOpen)
{

    EXPECT_EQ(stat_init, TNG_SUCCESS);
}


TEST_F(WaterVelsForcesTest, BoxShapeRead)
{
    int64_t             stride_length;
    tng_function_status read_stat;
    read_stat = tng_util_box_shape_read(traj, &box_shape, &stride_length);
    EXPECT_EQ(read_stat, TNG_SUCCESS);
    EXPECT_EQ(stride_length, 5000);
}

TEST_F(WaterVelsForcesTest, BoxShapeValues)
{
    int64_t             stride_length, i;
    const double        errtol = 1e-5;
    tng_function_status read_stat;
    read_stat = tng_util_box_shape_read(traj, &box_shape, &stride_length);
    EXPECT_EQ(read_stat, TNG_SUCCESS);
    EXPECT_EQ(stride_length, 5000);
    const std::vector<float> frame_0 = {
        2.87951e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 2.87951e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 2.87951e+00,
    };

    ASSERT_EQ(frame_0.size(), 9);
    for (i = 0; i < 9; i++)
    {
        EXPECT_NEAR(box_shape[i], frame_0[i], errtol);
    }
    const std::vector<float> frame_100 = { 2.89497e+00, 0.00000e+00, 0.00000e+00,
                                           0.00000e+00, 2.89497e+00, 0.00000e+00,
                                           0.00000e+00, 0.00000e+00, 2.89497e+00 };

    ASSERT_EQ(frame_100.size(), 9);
    for (i = 0; i < 9; i++)
    {
        EXPECT_NEAR(box_shape[909 - 9 + i], frame_100[i], errtol);
    }
}

TEST_F(WaterVelsForcesTest, PositionRead)
{
    int64_t             stride_length;
    tng_function_status read_stat;
    read_stat = tng_util_pos_read(traj, &positions, &stride_length);
    EXPECT_EQ(read_stat, TNG_SUCCESS);
    EXPECT_EQ(stride_length, 5000);
}

TEST_F(WaterVelsForcesTest, PositionValues)
{
    int64_t             start, end, stride_length, i;
    const double        errtol = 1e-5;
    tng_function_status read_stat;
    read_stat = tng_util_pos_read(traj, &positions, &stride_length);
    EXPECT_EQ(read_stat, TNG_SUCCESS);
    // xyz first 10 atoms frame 0
    // gmx dump frame 0
    // clang-format off
    const std::vector<float> frame_0_first_10_values  = {
    2.52700e+00,  2.61101e+00,  2.45398e+00,  
    2.50319e+00,  2.59390e+00,  2.54510e+00,  
    2.61687e+00,  2.57898e+00,  2.44623e+00, 
    1.09097e+00,  1.27301e+00,  1.99202e+00,   
    1.01457e+00,  1.23310e+00,  2.03366e+00,   
    1.13694e+00,  1.19976e+00,  1.95100e+00,   
    2.20399e+00,  1.37297e+00,  8.83017e-01,   
    2.13535e+00,  1.38523e+00,  9.48592e-01,   
    2.21780e+00,  1.46022e+00,  8.46139e-01,   
    1.10605e+00,  2.11799e+00,  5.61040e-01   
    };
    // clang-format on
    ASSERT_EQ(frame_0_first_10_values.size(), 30);
    // compare element by element
    for (i = 0; i < 30; i++)
    {
        EXPECT_NEAR(positions[i], frame_0_first_10_values[i], errtol);
    }
    // xyz last 10 atoms frame 100
    // gmx dump frame 100
    // clang-format off
    const std::vector<float> frame_100_last_10_values  = {
    7.98970e-01,  2.15481e+00,  2.75854e+00,
    6.32804e-01,  6.59262e-01,  1.12701e+00,
    5.47739e-01,  6.89158e-01,  1.09488e+00,
    6.16521e-01,  5.70554e-01,  1.15907e+00,
    5.33961e-01,  2.20212e+00,  6.22357e-02,
    4.79836e-01,  2.17921e+00,  1.37788e-01,
    4.79169e-01,  2.18181e+00,  2.88140e+00,
    5.76261e-01,  1.85258e+00,  1.69974e+00,
    6.60233e-01,  1.87443e+00,  1.74016e+00,
    5.79366e-01,  1.75766e+00,  1.68776e+00,

    };
    // clang-format on
    ASSERT_EQ(frame_100_last_10_values.size(), 30);
    // compare element by element
    for (i = 0; i < 30; i++)
    {
        EXPECT_NEAR(positions[818100 - 30 + i], frame_100_last_10_values[i], errtol);
    }
}


TEST_F(WaterVelsForcesTest, ForceRead)
{
    int64_t             stride_length;
    tng_function_status read_stat;
    read_stat = tng_util_force_read(traj, &forces, &stride_length);
    EXPECT_EQ(read_stat, TNG_SUCCESS);
}

TEST_F(WaterVelsForcesTest, ForceValues)
{
    int64_t             start, end, stride_length, i;
    const double        errtol = 1e-2; // larger as forces are larger
    tng_function_status read_stat;
    read_stat = tng_util_force_read(traj, &forces, &stride_length);
    EXPECT_EQ(read_stat, TNG_SUCCESS);
    // forces first 10 atoms frame 0
    // clang-format off
     const std::vector<float> frame_0_first_10_values  = {
    -4.35261e+02,  3.36017e+02, -9.38570e+02,
    -1.75984e+01, -2.44064e+02,  1.25406e+03,
     6.57882e+02, -2.07715e+02,  2.72886e+02,
     1.75474e+01,  1.57273e+03,  2.80544e+01,
    -5.30602e+02, -8.79351e+02,  2.76766e+02,
     7.45154e+01, -5.15662e+02, -3.61260e+02,
     4.70405e+02, -1.26065e+03, -2.68651e+02,
    -5.15954e+02,  5.19739e+02,  2.85984e+02,
    -3.90010e+02,  4.82308e+02,  2.96046e+00,
     1.23199e+03, -7.51883e+02, -6.58181e+02,
    };
    // clang-format on
    ASSERT_EQ(frame_0_first_10_values.size(), 30);
    // compare element by element
    for (i = 0; i < 30; i++)
    {
        EXPECT_NEAR(forces[i], frame_0_first_10_values[i], errtol);
    }
    // forces last 10 atoms frame 100
    // clang-format off
    const std::vector<float> frame_100_last_10_values  = {
    -4.49360e+02, -5.46652e+02,  5.24477e+02,
     1.27648e+03,  8.27699e+02,  2.98916e+01,
    -9.49143e+02, -3.13201e+02, -3.78830e+02,
    -5.04814e+02, -5.57331e+02, -6.48604e+01,
     1.24046e+03,  1.05411e+03,  4.06005e+02,
    -3.61442e+02, -5.29395e+02,  1.26982e+02,
    -4.76165e+02, -5.24370e+02, -3.48132e+02,
    -7.41153e+02,  1.19924e+01, -7.19316e+02,
     5.67011e+02,  6.64948e+01,  2.13465e+02,
     2.43871e+02, -4.09309e+02,  4.87609e+01,
    };
    // clang-format on
    ASSERT_EQ(frame_100_last_10_values.size(), 30);
    // compare element by element
    for (i = 0; i < 30; i++)
    {
        EXPECT_NEAR(forces[818100 - 30 + i], frame_100_last_10_values[i], errtol);
    }
}

TEST_F(WaterVelsForcesTest, VelRead)
{
    int64_t             stride_length;
    tng_function_status read_stat;
    read_stat = tng_util_vel_read(traj, &velocities, &stride_length);
    EXPECT_EQ(read_stat, TNG_SUCCESS);
}

TEST_F(WaterVelsForcesTest, VelValues)
{
    int64_t             start, end, stride_length, i;
    const double        errtol = 1e-5;
    tng_function_status read_stat;
    read_stat = tng_util_vel_read(traj, &velocities, &stride_length);
    EXPECT_EQ(read_stat, TNG_SUCCESS);
    // vels first 10 atoms frame 0
    // clang-format off
    const std::vector<float> frame_0_first_10_values  = {
     3.51496e-01,  7.29674e-01, -5.33343e-02,
     5.97873e-02, -1.00359e+00, -4.19582e-01,
     2.56209e-01,  5.52850e-01, -4.53435e-01,
    -1.09184e-02,  3.66412e-01, -4.85018e-01,
     9.26847e-01, -6.03737e-01,  3.67032e-01,
    -9.85010e-02,  1.09447e+00, -1.94833e+00,
    -4.60571e-02,  3.64507e-01, -2.01200e-01,
    -1.23912e+00, -3.46699e-01, -1.27041e+00,
     6.12738e-01,  7.64292e-01,  9.39986e-01,
    -6.34257e-02, -3.96772e-02, -4.55601e-01,
    };
    // clang-format on
    ASSERT_EQ(frame_0_first_10_values.size(), 30);
    // compare element by element
    for (i = 0; i < 30; i++)
    {
        EXPECT_NEAR(velocities[i], frame_0_first_10_values[i], errtol);
    }
    // vels last 10 atoms frame 100
    // clang-format off
    const std::vector<float> frame_100_last_10_values  = {
    -1.29712e+00,  1.89736e-01, -4.58020e-01,
    -2.24550e-01,  1.98991e-01, -7.18228e-01,
     9.92350e-02,  1.55654e-01, -1.64584e+00,
    -6.58128e-01,  4.26997e-01, -2.94439e-01,
    -2.47945e-01, -4.03298e-01,  2.42530e-01,
     3.88940e-01,  2.55276e-01,  9.15576e-01,
    -1.57709e+00,  5.61387e-01,  9.03308e-01,
    -5.50578e-01, -3.38237e-01, -9.82961e-02,
     4.52938e-01, -7.97070e-01, -1.83071e+00,
    -7.36810e-01, -2.02619e-01, -1.35719e+00,
    };
    // clang-format on
    ASSERT_EQ(frame_100_last_10_values.size(), 30);
    // compare element by element
    for (i = 0; i < 30; i++)
    {
        EXPECT_NEAR(velocities[818100 - 30 + i], frame_100_last_10_values[i], errtol);
    }
}