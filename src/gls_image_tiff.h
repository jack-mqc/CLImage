//
//  gls_image_tiff.hpp
//  CLImageTest
//
//  Created by Fabio Riccardi on 3/15/22.
//

#ifndef gls_image_tiff_hpp
#define gls_image_tiff_hpp

#include <functional>
#include <string>
#include <span>

namespace gls {

void read_tiff_file(const std::string& filename, int pixel_channels, int pixel_bit_depth,
                    std::function<bool(int width, int height)> image_allocator,
                    std::function<void(int tiff_bitspersample, int tiff_samplesperpixel, int row, int strip_height,
                                       uint8_t *tiff_buffer)> process_tiff_strip);

template <typename T>
bool write_tiff_file(const std::string& filename, int width, int height, int pixel_channels, int pixel_bit_depth,
                     /* int compression, */ std::function<T*(int row)> row_pointer);

}  // namespace gls

#endif /* gls_image_tiff_hpp */
