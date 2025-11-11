//
// Created by Harieasswar Lakshmidevi on 8/30/25.
//

#ifndef COMPRESS_TAGS_H
#define COMPRESS_TAGS_H

#include <string>

/**
 * Compresses a tags file using the compressed serialization format
 * @param input_tags_file Path to the input tags file
 * @param output_prefix Prefix for output files (will create _compressed.tags, _encoded_starts.bin, _bwt_intervals.bin)
 */
void compress_tags_file(const std::string& input_tags_file, const std::string& output_prefix);

#endif // COMPRESS_TAGS_H