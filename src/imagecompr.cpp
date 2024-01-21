#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION

#include <iostream>
#include <chrono>
#include <omp.h>

#include "include/stb_image.h"
#include "include/stb_image_write.h"

#include "utils.hpp"
#include "applications.hpp"
#include "fullMatrix.hpp"
#include "svd.hpp"

/**
 * Brief Main function to demonstrate image compression.
 * 
 * This program loads an input image in PNG format, compresses it using an SVD-based
 * method, and saves the compressed image as a PNG file.
 * 
 * return 0 if the program runs successfully, -1 on failure.
 */
int main()
{

    /**
     * Upload the input image in png format and
     * set the output file image
     */
    const char *filename = "input_image/luna_rossa.png";
    const char *outputFilename = "output_image/luna_rossa_out.png";

    /**
     * Create the object of the compression image application
    */
    APPLICATIONS obj;

    int Width, Height, channels;

    /**
     * Load the input image
    */
    stbi_uc *Image = stbi_load(filename, &Width, &Height, &channels, 0);


    if (Image == nullptr)
    {
        std::cerr << "Failed to load image: " << filename << std::endl;
        return -1;
    }

    /**
     * Create a matrix for compressed image data
    */
    Matrix Compressed(Height, Width);

    /**
     *  Perform image compression using SVD
    */

    Compressed = obj.image_compression(Image, channels, Height, Width);

    /**
     * Perform backward conversion to obtain the decompressed image
    */
    obj.backward_conversion(Image, Compressed, channels, Height, Width);

    /**
     * Save the decompressed image
    */
    stbi_write_png(outputFilename, Width, Height, channels, Image, Width * channels);

    /**
     * Free allocated memory
    */
    stbi_image_free(Image);


    return 0;
}