# ScanlineRenderer
Scanline renderer implemented using C++.
Last Updated: 2/2/2020

# Dependencies
 * Eigen: I used the matrix and vector classes from the Eigen library (http://eigen.tuxfamily.org/).

# Description
This is a scanline renderer, implemented using C++.

This code includes:

* Application of the viewport matrix
* Triangle rasterization
* Application of the orthographic projection matrix
* Application of the camera matrix
* Application of the perspective projection matrix
* Perspective divide
* Z-buffering
* Color interpolation

# Running the Code
The user can define values for the `eye`, `lookat`, and `up` vectors. An image named 'custom.ppm' will be produced using these values. 
Usage:

```
  make
  ./run [eye value 1] [eye value 2] [eye value 3] [lookat value 1] [lookat value 2] [lookat value 3] [up value 1] [up value 2] [up value 3]
```
An example output with `eye` = (1,1,1), `lookat` = (0.2,0.4,0), and `up` =(0,0,1) is included in the `examples` folder.

```
make
./run 1 1 1 0.2 0.4 0 0 0 1
```

Running the code with no command line arguments will produce 8 images demonstrating the result of each of the steps in the description above with `eye` = (1,1,1), `lookat` = (0,0,0), and `up` =(0,1,0).

```
  make
  ./run
````
# Potential Improvements/Additions
* Include texture mapping
* The transform() function applies a matrix to a set of vertices. In order to apply more than one, I have
to call the function consecutively with each transform. I would like to be able to multiply the matrices 
before calling the function, passing in the product and therefore only having to call the function once.
Unfortunately this results in the right picture but in the wrong position, and some other weird changes to the image,
so for now I just call it multiple times.
