// Copyright (c) 2020 Robert Vaser

#ifndef BIOSOUP_NUCLEIC_ACID_HPP_
#define BIOSOUP_NUCLEIC_ACID_HPP_

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <numeric>
#include <string>
#include <stdexcept>
#include <vector>
#include <iostream>

namespace biosoup {

constexpr static std::uint8_t kNucleotideCoder[] = {
    255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255,   0, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255,
    255,   0,   1 ,  1,   0, 255, 255,   2,
      3, 255, 255,   2, 255,   1,   0, 255,
    255, 255,   0,   1,   3,   3,   2,   0,
    255,   3, 255, 255, 255, 255, 255, 255,
    255,   0,   1,   1,   0, 255, 255,   2,
      3, 255, 255,   2, 255,   1,   0, 255,
    255, 255,   0,   1,   3,   3,   2,   0,
    255,   3, 255, 255, 255, 255, 255, 255
};

constexpr static char kNucleotideDecoder[] = {
    'A', 'C', 'G', 'T'
};

class NucleicAcid {
 public:
  NucleicAcid() = default;

  NucleicAcid(
      const std::string& name,
      const std::string& data)
      : NucleicAcid(
          name.c_str(), name.size(),
          data.c_str(), data.size()) {}

  NucleicAcid(
      const char* name, std::uint32_t name_len,
      const char* data, std::uint32_t data_len)
      : id(num_objects++),
        name(name, name_len),
        deflated_data(),
        block_quality(),
        inflated_len(data_len),
        is_reverse_complement(0) {
    deflated_data.reserve(data_len / 32. + .999);
    std::uint64_t block = 0;
    for (std::uint32_t i = 0; i < data_len; ++i) {
      std::uint64_t c = kNucleotideCoder[static_cast<std::uint8_t>(data[i])];
      if (c == 255ULL) {
        throw std::invalid_argument(
            "[biosoup::NucleicAcid::NucleicAcid] error: not a nucleotide");
      }
      block |= c << ((i << 1) & 63);
      if (((i + 1) & 31) == 0 || i == data_len - 1) {
        deflated_data.emplace_back(block);
        block = 0;
      }
    }
  }

  NucleicAcid(
      const std::string& name,
      const std::string& data,
      const std::string& quality)
      : NucleicAcid(
          name.c_str(), name.size(),
          data.c_str(), data.size(),
          quality.c_str(), quality.size()) {}

  NucleicAcid(
      const char* name, std::uint32_t name_len,
      const char* data, std::uint32_t data_len,
      const char* quality, std::uint32_t quality_len)
      : NucleicAcid(
          name, name_len,
          data, data_len) {
    
    deflated_quality.reserve(quality_len / 32. + .999);
    quality_levels.reserve(quality_len / 512. + .999);
    for (std::uint32_t i = 0; i < quality_len; i += 512) {         // for each window of 512 find the corresponding encoding and apply it
      std::uint32_t index_limit = std::min(i + 512, quality_len);  // upper bound of the window
      std::uint64_t quality_sum = 0;
      std::vector<std::int32_t> quality_freq(100, 0);              // counting sort of quality values
      for (std::uint32_t j = i; j < index_limit; ++j) {
        std::uint8_t curr_quality = quality[j] - '!';
        quality_freq[curr_quality]++;
        quality_sum += curr_quality;
      }
      std::vector<uint8_t> curr_levels;                            // current encoding values
      curr_levels.reserve(4);
      DecideQualityLevels(quality_freq, curr_levels, quality_sum, index_limit, i);
      //Encoding quality values of current window
      std::uint64_t block = 0;
      std::vector<std::int32_t> diffs;                             // vector of distances to each encoding
      diffs.reserve(4);
      for (std::uint32_t j = i; j < index_limit; ++j) {
        diffs.clear();
        std::uint8_t curr_quality = quality[j] - '!';
        for (auto x : curr_levels) diffs.emplace_back(std::abs(x - curr_quality));      // calculate distance from each encoding of current quality score 
        std::uint64_t c = std::min_element(diffs.begin(), diffs.end()) - diffs.begin(); // use the encoding closest to actual quality score
        block |= c << ((j << 1) & 63);                                                  // append current encoding to block
        if ( j == quality_len - 1 || ((j + 1) & 31) == 0 ) {                            // if block is full append it to deflated_quality
          deflated_quality.emplace_back(block);
          block = 0;
        }
      }
    }
  }
  

  NucleicAcid(const NucleicAcid&) = default;
  NucleicAcid& operator=(const NucleicAcid&) = default;

  NucleicAcid(NucleicAcid&&) = default;
  NucleicAcid& operator=(NucleicAcid&&) = default;

  ~NucleicAcid() = default;

  void DecideQualityLevels(std::vector<std::int32_t> &quality_freq,  
                           std::vector<uint8_t> &curr_levels,
                           std::uint64_t quality_sum,
                           std::uint32_t index_limit,
                           std::uint32_t i) {
    //Gathering information about the distribution shape on current window
    std::uint8_t min_q = std::find_if(
                            quality_freq.begin(), 
                            quality_freq.end(), 
                            [](std::int32_t x) -> bool { return x != 0; }
                            ) - quality_freq.begin();
    std::uint8_t max_q = quality_freq.rend() - std::find_if(
                                                  quality_freq.rbegin(), 
                                                  quality_freq.rend(), 
                                                  [](std::int32_t x) -> bool { return x != 0; }
                                                  ) - 1;
    std::uint8_t avg_q = quality_sum / (static_cast<std::int32_t>(index_limit) - static_cast<std::int32_t>(i));
    std::uint8_t mod_q = std::max_element(
                            quality_freq.begin(), 
                            quality_freq.end(), 
                            [](std::int32_t a, std::int32_t b) -> bool { return a < b; }
                            ) - quality_freq.begin();
    //Deciding quality discretization levels for current window
    double quarter = (max_q - min_q) / 4.0;
    if (quarter < 1.0) quarter = 1.0;
    std::int32_t num_of_upper_levels;
    std::int32_t num_of_lower_levels;
    if (mod_q > avg_q) {              // if distribution is negatively skewed
      num_of_upper_levels = static_cast<std::int32_t>((max_q - mod_q) / quarter);
      num_of_lower_levels = 4 - num_of_upper_levels - 1;
    } else if (mod_q < avg_q) {      // if distribution is positively skewed
      num_of_lower_levels = static_cast<std::int32_t>((mod_q - min_q) / quarter);
      num_of_upper_levels = 4 - num_of_lower_levels - 1;
    } else {                         // if distribution is not skewed
      num_of_lower_levels = 1;       // the mode value will be considered a lower level
      num_of_upper_levels = 2;
    }
    std::int32_t upper_step = (max_q - mod_q) / (num_of_upper_levels + 1);
    std::int32_t lower_step = (mod_q - min_q) / (num_of_lower_levels + 1);
    std::uint32_t discretization_levels = 0;  // to be appended to quality levels
    std::uint32_t curr_level = min_q;
    auto process_level = [&curr_levels, &discretization_levels] (std::uint32_t level) {
      curr_levels.emplace_back(static_cast<uint8_t>(level));
      discretization_levels <<= 8;
      discretization_levels |= level;
    };
    for (int32_t j = 0; j < num_of_lower_levels; j++) {
      curr_level += lower_step;
      process_level(curr_level);
    }
    curr_level = mod_q;
    process_level(curr_level);
    for (int32_t j = 0; j < num_of_upper_levels; j++) {
      curr_level += upper_step;
      process_level(curr_level);
    }
    quality_levels.push_back(discretization_levels);       // append encoding for current window to quality_levels vector
    std::reverse(curr_levels.begin(), curr_levels.end());  // index of encoding in curr_levels must be proportional to 
  }                                                        // how much discretization_levels have to be shifted to the right to get the same encoding

  std::uint64_t Code(std::uint32_t i) const {
    std::uint64_t x = 0;
    if (is_reverse_complement) {
      i = inflated_len - i - 1;
      x = 3;
    }
    return ((deflated_data[i >> 5] >> ((i << 1) & 63)) & 3) ^ x;
  }

  std::uint8_t Score(std::uint32_t i) const {
    if (is_reverse_complement) {
      i = inflated_len - i - 1;
    }
    std::int32_t index = (deflated_quality[i >> 5] >> ((i << 1) & 63)) & 3; // offset of the decoding
    return static_cast<std::uint8_t>((quality_levels[i >> 9] >> (index * 8)) & 255);
  }

  std::string InflateData(std::uint32_t i = 0, std::uint32_t len = -1) const {
    if (i >= inflated_len) {
      return std::string{};
    }
    len = std::min(len, inflated_len - i);

    std::string dst{};
    dst.reserve(len);
    for (; len; ++i, --len) {
      dst += kNucleotideDecoder[Code(i)];
    }
    return dst;
  }

  std::string InflateQuality(std::uint32_t i = 0, std::uint32_t len = -1) const {  // NOLINT
    if (deflated_quality.empty() || i >= inflated_len) {
      return std::string{};
    }
    len = std::min(len, inflated_len - i);

    std::string dst{};
    dst.reserve(len);
    for (; len; ++i, --len) {
      dst += Score(i) + '!';
    }
    return dst;
  }

  void ReverseAndComplement() {   // Watson-Crick base pairing
    is_reverse_complement ^= 1;
  }

  static std::atomic<std::uint32_t> num_objects;

  std::uint32_t id;  // (optional) initialize num_objects to 0
  std::string name;
  std::vector<std::uint64_t> deflated_data;
  std::vector<std::uint64_t> deflated_quality; // each Phred quality score encoded with 2 bits
  std::vector<std::uint8_t> block_quality;  // (optional) Phred quality scores
  std::vector<std::uint32_t> quality_levels; // 2 bit decodings for each window of 512 quality scores
  std::uint32_t inflated_len;
  bool is_reverse_complement;
};

}  // namespace biosoup

#endif  // BIOSOUP_NUCLEIC_ACID_HPP_
