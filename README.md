# ScanlineRenderer
Scanline renderer implemented using C++.
Last Updated: 2/2/2020

# Dependencies
 * Eigen: I used the matrix and vector classes from the Eigen library (http://eigen.tuxfamily.org/).

# Description
This is a scanline renderer, implemented using C++.

This code includes:

* Triangle rasterization
* Application of the viewport matrix
* Application of the orthographic projection matrix
* Application of the perspective projection matrix
* Perspective divide
* Z-buffering
* Color interpolation

# Potential Improvements/Additions
* Include texture mapping
* The transform() function applies a matrix to a set of vertices. In order to apply more than one, I have
to call the function consecutively with each transform. I would like to be able to multiply the matrices 
before calling the function, passing in the product and therefore only having to call the function once.
Unfortunately this results in the right picture but in the wrong position, and some other weird changes to the image,
so for now I just call it multiple times.
