# openvdb_catkin
A catkin wrapper for http://www.openvdb.org/

Tested on Ubuntu 18.04, GCC 7.4.0

Includes the catkin package `openvdb_interface_test` that links against the
catkinized version and runs several of the tutorial examples to check that
everything is in order.

#### System Dependencies

```bash
sudo apt install -y libboost-iostreams-dev \
                    libboost-system-dev \
                    libtbb-dev \
                    libilmbase-dev \
                    libopenexr-dev
```
