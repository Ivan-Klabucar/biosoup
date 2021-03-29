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
    std::uint64_t quality_sum = 0;
    std::uint8_t min_quality = quality[0];
    std::uint8_t max_quality = quality[0];
    for (std::uint32_t i = 0; i < data_len; ++i) {
      std::uint8_t curr_quality = quality[i] - '!';
      quality_sum += curr_quality;
      if (curr_quality > max_quality) max_quality = curr_quality;
      if (curr_quality < min_quality) min_quality = curr_quality;
    }
    average_quality = quality_sum / quality_len;
    DecideQualityLevels(min_quality, max_quality);
  }
  

  NucleicAcid(const NucleicAcid&) = default;
  NucleicAcid& operator=(const NucleicAcid&) = default;

  NucleicAcid(NucleicAcid&&) = default;
  NucleicAcid& operator=(NucleicAcid&&) = default;

  ~NucleicAcid() = default;

  void DecideQualityLevels(std::uint8_t min_quality, std::uint8_t max_quality) {
    std::int32_t range = max_quality - min_quality;
    std::int32_t quarter = range / 4;
    std::int32_t num_of_upper_levels = (max_quality - average_quality) / quarter;
    if (num_of_upper_levels == 4) num_of_upper_levels = 3;
    std::int32_t num_of_lower_levels = 4 - num_of_upper_levels;
    std::int32_t upper_step = (max_quality - average_quality) / (num_of_upper_levels + 1);
    std::int32_t lower_step = (average_quality - min_quality) / num_of_lower_levels;
    std::int8_t curr_level = min_quality;
    while(num_of_lower_levels > 0) {
      curr_level += lower_step; 
      quality_levels.push_back(curr_level);
      num_of_lower_levels--;
    }
    curr_level = average_quality;
    while(num_of_upper_levels > 0) {
      curr_level += upper_step; 
      quality_levels.push_back(curr_level);
      num_of_upper_levels--;
    }
    std::cout << "Quality levels: " << std::endl;
    for(auto x : quality_levels) {
      std::cout << x << std::endl;
    }
  }

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
    return block_quality[i >> 6];
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
    if (block_quality.empty() || i >= inflated_len) {
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
  std::vector<std::uint64_t> deflated_quality;
  std::vector<std::uint8_t> block_quality;  // (optional) Phred quality scores
  std::vector<std::int8_t> quality_levels;
  std::uint8_t average_quality;
  std::uint32_t inflated_len;
  bool is_reverse_complement;
};

}  // namespace biosoup

#endif  // BIOSOUP_NUCLEIC_ACID_HPP_
