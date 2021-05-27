//
// Created by soliareofastora on 23.03.2021.
//

#ifndef LINESOCCER__BIT_SET_HPP_
#define LINESOCCER__BIT_SET_HPP_

struct BitSet {
  unsigned char *data;
  unsigned int size;

  BitSet(unsigned int N) {
    size = N / 8;
    if (N % 8 > 0) {
      size += 1;
    }
    data = new unsigned char[size];
  }

  ~BitSet() {
    delete[] data;
  }

  void set_bit(unsigned int i, bool value) const {
    if (value)
      data[i / 8] |= true << (i % 8);
    else
      data[i / 8] &= ~(true << (i % 8));
  }

  bool get_bit(unsigned int i) const {
    return (data[i / 8] >> (i % 8)) & true;
  }

  void set_char(unsigned int i, unsigned char value) const {
    data[i] = value;
  }
};

#endif  // LINESOCCER__BIT_SET_HPP_
