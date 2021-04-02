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
#include <map>
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
    std::map<uint8_t, int32_t> quality_freq;
    for (std::uint32_t i = 0; i < quality_len; ++i) {
      std::uint8_t curr_quality = quality[i] - '!';
      quality_freq[curr_quality]++;
      quality_sum += curr_quality;
    }
    std::uint8_t min_q = quality_freq.begin()->first;
    std::uint8_t max_q = quality_freq.rbegin()->first;
    std::uint8_t mod_q = std::max_element(
        quality_freq.begin(), 
        quality_freq.end(), 
        [] (const std::pair<uint8_t, int32_t>& a, const std::pair<uint8_t, int32_t>& b)-> bool { 
          return a.second < b.second; 
        })->first;
    std::uint8_t avg_q = quality_sum / quality_len;
    DecideQualityLevels(min_q, max_q, avg_q, mod_q);
    deflated_quality.reserve(quality_len / 32. + .999);
    std::uint64_t block = 0;
    for (std::uint32_t i = 0; i < quality_len; ++i) {
      std::uint8_t curr_quality = quality[i] - '!';
      std::vector<std::int32_t> diffs;
      std::for_each(
          quality_levels.begin(), 
          quality_levels.end(),
          [&diffs, &curr_quality] (std::uint8_t &x) { diffs.push_back(std::abs(x - curr_quality)); });
      std::uint64_t c = std::min_element(diffs.begin(), diffs.end()) - diffs.begin();
      block |= c << ((i << 1) & 63);
      if (((i + 1) & 31) == 0 || i == data_len - 1) {
        deflated_quality.emplace_back(block);
        block = 0;
      }
    }
  }
  

  NucleicAcid(const NucleicAcid&) = default;
  NucleicAcid& operator=(const NucleicAcid&) = default;

  NucleicAcid(NucleicAcid&&) = default;
  NucleicAcid& operator=(NucleicAcid&&) = default;

  ~NucleicAcid() = default;

  void DecideQualityLevels(std::uint8_t min_q, std::uint8_t max_q, std::uint8_t avg_q, std::uint8_t mod_q) {
    std::int32_t range = max_q - min_q;
    std::int32_t quarter = range / 4;
    if( quarter == 0 ) quarter = 1;
    std::int32_t num_of_upper_levels = (max_q - avg_q) / quarter;
    if ( num_of_upper_levels == 4 ) num_of_upper_levels = 3;
    std::int32_t num_of_lower_levels = 4 - num_of_upper_levels;
    std::int32_t upper_step = (max_q - avg_q) / (num_of_upper_levels + 1);
    std::int32_t lower_step = (avg_q - min_q) / num_of_lower_levels;
    std::int32_t curr_level = min_q;
    // std::cout << "min q: "<< (std::int32_t) min_q << std::endl;
    // std::cout << "max q: "<< (std::int32_t) max_q << std::endl;
    // std::cout << "upper s: "<< (std::int32_t) upper_step << std::endl;
    // std::cout << "lowwer s: "<< (std::int32_t) lower_step << std::endl;
    // std::cout << "avg: "<< (std::int32_t) average_quality << std::endl;
    while(num_of_lower_levels > 0) {
      curr_level += lower_step; 
      quality_levels.push_back(curr_level);
      num_of_lower_levels--;
    }
    curr_level = avg_q;
    while(num_of_upper_levels > 0) {
      curr_level += upper_step; 
      quality_levels.push_back(curr_level);
      num_of_upper_levels--;
    }
    // std::cout << "Quality levels: " << std::endl;
    // for(auto x : quality_levels) {
    //   std::cout << (std::int32_t)x << std::endl;
    // }
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
    return quality_levels[(deflated_quality[i >> 5] >> ((i << 1) & 63)) & 3];
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
  std::vector<std::uint64_t> deflated_quality;
  std::vector<std::uint8_t> block_quality;  // (optional) Phred quality scores
  std::vector<std::uint8_t> quality_levels;
  std::uint32_t inflated_len;
  bool is_reverse_complement;
};

}  // namespace biosoup

#endif  // BIOSOUP_NUCLEIC_ACID_HPP_
