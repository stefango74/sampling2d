sample2d reads a SVG vector file containing cubic Bezier curves, samples the curve, and meshes it.
Outputs are a .PLY curve, a .STL mesh, and .SVG vector files for both as well as some statistics on the command line

Build:
install CGAL, eigen and ANN libraries
On Linux this can be done by:
sudo apt install libcgal-dev libeigen3-dev libann-dev
make

Usage:
sample2d input.svg outfile-basename epsilon [mindist] [maxdist] [maxedgelen]
The optional arguments are in function of the bounding box diagonal of the control points

There are some sample .SVG files in the subdirectory ./data
I.e., you can try 'sample2d data/simple.svg test 0.33' to generate a 0.33-sampling of a simple curve
