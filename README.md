ctvr: Contour Tree Volume Renderer

This is minimally-maintained code that I've been dragging around with me for years. I'm not sure if it builds anymore. It's been years since I've had Qt installed. But Tourtre.hpp and Trilinear.hpp are very solid nuggets of code for computing contour trees (the former) of 3d grayscale images (the latter). The rest of the code is a pretty good example for how to use those two libraries. If you are thinking of using libtourtre (http://github.com/sedillard/libtourtre) but are put off by the minimalist ANSI C style, and want something more C++-templately, then try Tourtre.hpp. See ContourTree.cpp for example usage.

The method is described in the following publication:

Weber, G. H., Dillard, S. E., Carr, H., Pascucci, V., & Hamann, B. (2007). Topology-controlled volume rendering. Visualization and Computer Graphics, IEEE Transactions on, 13(2), 330-341.
http://escholarship.org/uc/item/6bt8m788

