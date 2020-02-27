#include <fstream>
#include <iostream>
#include <string>

#include <glog/logging.h>
#include <gtest/gtest.h>
#include <openvdb/io/Stream.h>
#include <openvdb/openvdb.h>
#include <openvdb/tools/ChangeBackground.h>
#include <openvdb/tools/LevelSetSphere.h>

class OpenVDBInterfaceTest : public ::testing::Test {
public:
protected:
  virtual void SetUp() {}

  virtual void TearDown() {}

public:
};

// Source:
// https://www.openvdb.org/documentation/doxygen/codeExamples.html#sPopulatingGrids
/////////////////////////////////////////////////////////////////////////////////////////
// Populate the given grid with a narrow-band level set representation of a
// sphere. The width of the narrow band is determined by the grid's background
// value. (Example code only; use tools::createSphereSDF() in production.)
template <class GridType>
void makeSphere(GridType &grid, float radius, const openvdb::Vec3f &c) {
  using ValueT = typename GridType::ValueType;
  // Distance value for the constant region exterior to the narrow band
  const ValueT outside = grid.background();
  // Distance value for the constant region interior to the narrow band
  // (by convention, the signed distance is negative in the interior of
  // a level set)
  const ValueT inside = -outside;
  // Use the background value as the width in voxels of the narrow band.
  // (The narrow band is centered on the surface of the sphere, which
  // has distance 0.)
  int padding = int(openvdb::math::RoundUp(openvdb::math::Abs(outside)));
  // The bounding box of the narrow band is 2*dim voxels on a side.
  int dim = int(radius + padding);
  // Get a voxel accessor.
  typename GridType::Accessor accessor = grid.getAccessor();
  // Compute the signed distance from the surface of the sphere of each
  // voxel within the bounding box and insert the value into the grid
  // if it is smaller in magnitude than the background value.
  openvdb::Coord ijk;
  int &i = ijk[0], &j = ijk[1], &k = ijk[2];
  for (i = c[0] - dim; i < c[0] + dim; ++i) {
    const float x2 = openvdb::math::Pow2(i - c[0]);
    for (j = c[1] - dim; j < c[1] + dim; ++j) {
      const float x2y2 = openvdb::math::Pow2(j - c[1]) + x2;
      for (k = c[2] - dim; k < c[2] + dim; ++k) {
        // The distance from the sphere surface in voxels
        const float dist =
            openvdb::math::Sqrt(x2y2 + openvdb::math::Pow2(k - c[2])) - radius;
        // Convert the floating-point distance to the grid's value type.
        ValueT val = ValueT(dist);
        // Only insert distances that are smaller in magnitude than
        // the background value.
        if (val < inside || outside < val)
          continue;
        // Set the distance for voxel (i,j,k).
        accessor.setValue(ijk, val);
      }
    }
  }
  // Propagate the outside/inside sign information from the narrow band
  // throughout the grid.
  openvdb::tools::signedFloodFill(grid.tree());
}

// Source:
// https://www.openvdb.org/documentation/doxygen/codeExamples.html#sHelloWorld
///////////////////////////////////////////////////////////////////////////////////
TEST(OpenVDBInterfaceTest, HelloWorld) {
  // Initialize the OpenVDB library.  This must be called at least
  // once per program and may safely be called multiple times.
  openvdb::initialize();
  // Create an empty floating-point grid with background value 0.
  openvdb::FloatGrid::Ptr grid = openvdb::FloatGrid::create();
  LOG(INFO) << "Testing random access:" << std::endl;
  // Get an accessor for coordinate-based access to voxels.
  openvdb::FloatGrid::Accessor accessor = grid->getAccessor();
  // Define a coordinate with large signed indices.
  openvdb::Coord xyz(1000, -200000000, 30000000);
  // Set the voxel value at (1000, -200000000, 30000000) to 1.
  accessor.setValue(xyz, 1.0);
  // Verify that the voxel value at (1000, -200000000, 30000000) is 1.
  LOG(INFO) << "Grid" << xyz << " = " << accessor.getValue(xyz) << std::endl;
  // Reset the coordinates to those of a different voxel.
  xyz.reset(1000, 200000000, -30000000);
  // Verify that the voxel value at (1000, 200000000, -30000000) is
  // the background value, 0.
  LOG(INFO) << "Grid" << xyz << " = " << accessor.getValue(xyz) << std::endl;
  // Set the voxel value at (1000, 200000000, -30000000) to 2.
  accessor.setValue(xyz, 2.0);
  // Set the voxels at the two extremes of the available coordinate space.
  // For 32-bit signed coordinates these are (-2147483648, -2147483648,
  // -2147483648) and (2147483647, 2147483647, 2147483647).
  accessor.setValue(openvdb::Coord::min(), 3.0f);
  accessor.setValue(openvdb::Coord::max(), 4.0f);
  LOG(INFO) << "Testing sequential access:" << std::endl;
  // Print all active ("on") voxels by means of an iterator.
  for (openvdb::FloatGrid::ValueOnCIter iter = grid->cbeginValueOn(); iter;
       ++iter) {
    LOG(INFO) << "Grid" << iter.getCoord() << " = " << *iter << std::endl;
  }
}

// Source:
// https://www.openvdb.org/documentation/doxygen/codeExamples.html#sAllocatingGrids
///////////////////////////////////////////////////////////////////////////////////
TEST(OpenVDBInterfaceTest, AllocatingGrids) {
  openvdb::initialize();
  // Create a shared pointer to a newly-allocated grid of a built-in type:
  // in this case, a FloatGrid, which stores one single-precision floating point
  // value per voxel.  Other built-in grid types include BoolGrid, DoubleGrid,
  // Int32Grid and Vec3SGrid (see openvdb.h for the complete list).
  // The grid comprises a sparse tree representation of voxel data,
  // user-supplied metadata and a voxel space to world space transform,
  // which defaults to the identity transform.
  openvdb::FloatGrid::Ptr grid =
      openvdb::FloatGrid::create(/*background value=*/2.0);
  // Populate the grid with a sparse, narrow-band level set representation
  // of a sphere with radius 50 voxels, located at (1.5, 2, 3) in index space.
  makeSphere(*grid, /*radius=*/50.0, /*center=*/openvdb::Vec3f(1.5, 2, 3));
  // Associate some metadata with the grid.
  grid->insertMeta("radius", openvdb::FloatMetadata(50.0));
  // Associate a scaling transform with the grid that sets the voxel size
  // to 0.5 units in world space.
  grid->setTransform(
      openvdb::math::Transform::createLinearTransform(/*voxel size=*/0.5));
  // Identify the grid as a level set.
  grid->setGridClass(openvdb::GRID_LEVEL_SET);
  // Name the grid "LevelSetSphere".
  grid->setName("LevelSetSphere");
  // Create a VDB file object.
  openvdb::io::File file("mygrids.vdb");
  // Add the grid pointer to a container.
  openvdb::GridPtrVec grids;
  grids.push_back(grid);
  // Write out the contents of the container.
  file.write(grids);
  file.close();
}
TEST(OpenVDBInterfaceTest, AllocatingGridsShortVersion) {
  openvdb::initialize();
  // Create a FloatGrid and populate it with a narrow-band
  // signed distance field of a sphere.
  openvdb::FloatGrid::Ptr grid =
      openvdb::tools::createLevelSetSphere<openvdb::FloatGrid>(
          /*radius=*/50.0, /*center=*/openvdb::Vec3f(1.5, 2, 3),
          /*voxel size=*/0.5, /*width=*/4.0);
  // Associate some metadata with the grid.
  grid->insertMeta("radius", openvdb::FloatMetadata(50.0));
  // Name the grid "LevelSetSphere".
  grid->setName("LevelSetSphere");
  // Create a VDB file object and write out the grid.
  openvdb::io::File("mygrids.vdb").write({grid});
}

// Source:
// https://www.openvdb.org/documentation/doxygen/codeExamples.html#sPopulatingGrids
///////////////////////////////////////////////////////////////////////////////////
TEST(OpenVDBInterfaceTest, PopulatingGrids) {
  openvdb::initialize();
  // Create a VDB file object.
  openvdb::io::File file("mygrids.vdb");
  // Open the file.  This reads the file header, but not any grids.
  file.open();
  // Loop over all grids in the file and retrieve a shared pointer
  // to the one named "LevelSetSphere".  (This can also be done
  // more simply by calling file.readGrid("LevelSetSphere").)
  openvdb::GridBase::Ptr baseGrid;
  for (openvdb::io::File::NameIterator nameIter = file.beginName();
       nameIter != file.endName(); ++nameIter) {
    // Read in only the grid we are interested in.
    if (nameIter.gridName() == "LevelSetSphere") {
      baseGrid = file.readGrid(nameIter.gridName());
    } else {
      LOG(INFO) << "skipping grid " << nameIter.gridName() << std::endl;
    }
  }
  file.close();
  // From the example above, "LevelSetSphere" is known to be a FloatGrid,
  // so cast the generic grid pointer to a FloatGrid pointer.
  openvdb::FloatGrid::Ptr grid =
      openvdb::gridPtrCast<openvdb::FloatGrid>(baseGrid);
  // Convert the level set sphere to a narrow-band fog volume, in which
  // interior voxels have value 1, exterior voxels have value 0, and
  // narrow-band voxels have values varying linearly from 0 to 1.
  const float outside = grid->background();
  const float width = 2.0 * outside;
  // Visit and update all of the grid's active values, which correspond to
  // voxels on the narrow band.
  for (openvdb::FloatGrid::ValueOnIter iter = grid->beginValueOn(); iter;
       ++iter) {
    float dist = iter.getValue();
    iter.setValue((outside - dist) / width);
  }
  // Visit all of the grid's inactive tile and voxel values and update the
  // values that correspond to the interior region.
  for (openvdb::FloatGrid::ValueOffIter iter = grid->beginValueOff(); iter;
       ++iter) {
    if (iter.getValue() < 0.0) {
      iter.setValue(1.0);
      iter.setValueOff();
    }
  }
  // Set exterior voxels to 0.
  openvdb::tools::changeBackground(grid->tree(), 0.0);
}

// Source:
// https://www.openvdb.org/documentation/doxygen/codeExamples.html#sStreamIO
///////////////////////////////////////////////////////////////////////////////////
TEST(OpenVDBInterfaceTest, StreamIO) {
  openvdb::initialize();
  // Create a VDB file object.
  openvdb::io::File file("mygrids.vdb");
  // Open the file.  This reads the file header, but not any grids.
  file.open();

  openvdb::GridBase::Ptr baseGrid;
  for (openvdb::io::File::NameIterator nameIter = file.beginName();
       nameIter != file.endName(); ++nameIter) {
    // Read in only the grid we are interested in.
    if (nameIter.gridName() == "LevelSetSphere") {
      baseGrid = file.readGrid(nameIter.gridName());
    } else {
      LOG(INFO) << "skipping grid " << nameIter.gridName() << std::endl;
    }
  }
  file.close();

  openvdb::GridPtrVecPtr grids(new openvdb::GridPtrVec);
  grids->push_back(baseGrid);

  // Stream the grids to a string.
  std::ostringstream ostr(std::ios_base::binary);
  openvdb::io::Stream(ostr).write(*grids);
  // Stream the grids to a file.
  std::ofstream ofile("mygrids.vdb", std::ios_base::binary);
  openvdb::io::Stream(ofile).write(*grids);
  // Stream in grids from a string.
  // Note that io::Stream::getGrids() returns a shared pointer
  // to an openvdb::GridPtrVec.
  std::istringstream istr(ostr.str(), std::ios_base::binary);
  openvdb::io::Stream strm(istr);
  grids = strm.getGrids();
  // Stream in grids from a file.
  std::ifstream ifile("mygrids.vdb", std::ios_base::binary);
  grids = openvdb::io::Stream(ifile).getGrids();
}

int main(int argc, char **argv) {
  google::InitGoogleLogging(argv[0]);
  google::InstallFailureSignalHandler();
  testing::InitGoogleTest(&argc, argv);
  google::ParseCommandLineFlags(&argc, &argv, true);
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  FLAGS_alsologtostderr = true;
  FLAGS_colorlogtostderr = true;
  FLAGS_v = 1;
  return RUN_ALL_TESTS();
}
