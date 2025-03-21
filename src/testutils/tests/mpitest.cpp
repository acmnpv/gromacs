/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2016- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
/*! \internal \file
 * \brief
 * Tests for infrastructure for running tests under MPI.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_testutils
 */
#include "gmxpre.h"

#include "testutils/mpitest.h"

#include "config.h"

#include <string>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/utility/basenetwork.h"
#include "gromacs/utility/gmxmpi.h"

#include "testutils/refdata.h"

namespace gmx
{
namespace test
{
namespace
{

class MpiSelfTest : public ::testing::Test
{
public:
    //! Whether each rank participated, relevant only on rank 0
    std::vector<int> reached_;
};

TEST_F(MpiSelfTest, Runs)
{
    GMX_MPI_TEST(RequireMinimumRankCount<2>);
    if (gmx_node_rank() == 0)
    {
        reached_.resize(getNumberOfTestMpiRanks(), 0);
    }
    // Needed for thread-MPI so that we resize the buffer before we
    // fill it on non-main ranks.
    MPI_Barrier(MPI_COMM_WORLD);
    int value = 1;
    MPI_Gather(&value, 1, MPI_INT, reached_.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (gmx_node_rank() == 0)
    {
        EXPECT_THAT(reached_, testing::Each(value));
    }
}

class MpiRefDataTest : public ::testing::Test
{
public:
};

TEST_F(MpiRefDataTest, UsesDifferentReferenceDataOnEachRank)
{
    GMX_MPI_TEST(RequireRankCount<2>);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // Test that one can initialize test reference data as a different
    // rank.
    TestReferenceData    refData(rank);
    TestReferenceChecker rootChecker(refData.rootChecker());
    // Test that reference data can actually differ between ranks.
    rootChecker.checkInteger(rank, "MPI rank");
}

TEST_F(MpiRefDataTest, UsesDifferentReferenceDataOnEachRankWithCustomName)
{
    GMX_MPI_TEST(RequireRankCount<2>);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    TestReferenceData    refData("CustomName.xml", rank);
    TestReferenceChecker rootChecker(refData.rootChecker());
    rootChecker.checkInteger(rank, "MPI rank");
}

} // namespace
} // namespace test
} // namespace gmx
