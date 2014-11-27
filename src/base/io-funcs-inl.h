// base/io-funcs-inl.h

// Copyright 2009-2011  Microsoft Corporation;  Saarland University;
//                      Jan Silovsky;   Yanmin Qian

// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at

//  http://www.apache.org/licenses/LICENSE-2.0

// THIS CODE IS PROVIDED *AS IS* BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
// KIND, EITHER EXPRESS OR IMPLIED, INCLUDING WITHOUT LIMITATION ANY IMPLIED
// WARRANTIES OR CONDITIONS OF TITLE, FITNESS FOR A PARTICULAR PURPOSE,
// MERCHANTABLITY OR NON-INFRINGEMENT.
// See the Apache 2 License for the specific language governing permissions and
// limitations under the License.

#ifndef KALDI_BASE_IO_FUNCS_INL_H_
#define KALDI_BASE_IO_FUNCS_INL_H_ 1

// Do not include this file directly.  It is included by base/io-funcs.h

#include <limits>
#include <vector>

namespace kaldi {

// Template that covers integers.
template<class T>  void WriteBasicType(std::ostream &os,
                                       bool binary, T t) {
  // Compile time assertion that this is not called with a wrong type.
  KALDI_ASSERT_IS_INTEGER_TYPE(T);
  if (binary) {
    char len_c = (std::numeric_limits<T>::is_signed ? 1 :  -1)
        * static_cast<char>(sizeof(t));
    os.put(len_c);
    os.write(reinterpret_cast<const char *>(&t), sizeof(t));
  } else {
    if (sizeof(t) == 1)
      os << static_cast<int16>(t) << " ";
    else
      os << t << " ";
  }
  if (os.fail()) {
    throw std::runtime_error("Write failure in WriteBasicType.");
  }
}

// Template that covers integers.
template<class T> inline void ReadBasicType(std::istream &is,
                                            bool binary, T *t) {
#ifdef KALDI_PARANOID
  assert(t != NULL);
#endif
  // Compile time assertion that this is not called with a wrong type.
  KALDI_ASSERT_IS_INTEGER_TYPE(T);
  if (binary) {
    char len_c = is.get(), len_c_expected
        = (std::numeric_limits<T>::is_signed ? 1 :  -1)
        * static_cast<char>(sizeof(*t));

    if (len_c !=  len_c_expected) {
      KALDI_ERR << "ReadBasicType: did not get expected integer type, "
                << static_cast<int>(len_c)
                << " vs. " << static_cast<int>(len_c_expected)
                << ".  You can change this code to successfully"
                << " read it later, if needed.";
      // insert code here to read "wrong" type.  Might have a switch statement.
    }
    is.read(reinterpret_cast<char *>(t), sizeof(*t));
  } else {
    if (sizeof(*t) == 1) {
      int16 i;
      is >> i;
      *t = i;
    } else {
      is >> *t;
    }
  }
  if (is.fail()) {
    KALDI_ERR << "Read failure in ReadBasicType, file position is "
              << is.tellg() << ", next char is " << is.peek();
  }
}


template<class T> inline void WriteIntegerVector(std::ostream &os, bool binary,
                                                 const std::vector<T> &v) {
  // Compile time assertion that this is not called with a wrong type.
  KALDI_ASSERT_IS_INTEGER_TYPE(T);
  if (binary) {
    char sz = sizeof(T);  // this is currently just a check.
    os.write(&sz, 1);
    int32 vecsz = static_cast<int32>(v.size());
    assert((size_t)vecsz == v.size());
    os.write(reinterpret_cast<const char *>(&vecsz), sizeof(vecsz));
    if (vecsz != 0) {
      os.write(reinterpret_cast<const char *>(&(v[0])), sizeof(T)*vecsz);
    }
  } else {
    // focus here is on prettiness of text form rather than
    // efficiency of reading-in.
    // reading-in is dominated by low-level operations anyway:
    // for efficiency use binary.
    os << "[ ";
    typename std::vector<T>::const_iterator iter = v.begin(), end = v.end();
    for (; iter != end; ++iter) {
      if (sizeof(T) == 1)
        os << static_cast<int16>(*iter) << " ";
      else
        os << *iter << " ";
    }
    os << "]\n";
  }
  if (os.fail()) {
    throw std::runtime_error("Write failure in WriteIntegerType.");
  }
}


template<class T> inline void ReadIntegerVector(std::istream &is,
                                                bool binary,
                                                std::vector<T> *v) {
  KALDI_ASSERT_IS_INTEGER_TYPE(T);
  assert(v != NULL);
  if (binary) {
    int sz = is.peek();
    if (sz == sizeof(T)) {
      is.get();
    } else {  // this is currently just a check.
      KALDI_ERR << "ReadIntegerVector: expected to see type of size "
                << sizeof(T) << ", saw instead " << sz << ", at file position "
                << is.tellg();
    }
    int32 vecsz;
    is.read(reinterpret_cast<char *>(&vecsz), sizeof(vecsz));
    if (is.fail() || vecsz < 0) goto bad;
    v->resize(vecsz);
    if (vecsz > 0) {
      is.read(reinterpret_cast<char *>(&((*v)[0])), sizeof(T)*vecsz);
    }
  } else {
    std::vector<T> tmp_v;  // use temporary so v doesn't use extra memory
                           // due to resizing.
    is >> std::ws;
    if (is.peek() != static_cast<int>('[')) {
      KALDI_ERR << "ReadIntegerVector: expected to see [, saw "
                << is.peek() << ", at file position " << is.tellg();
    }
    is.get();  // consume the '['.
    is >> std::ws;  // consume whitespace.
    while (is.peek() != static_cast<int>(']')) {
      if (sizeof(T) == 1) {  // read/write chars as numbers.
        int16 next_t;
        is >> next_t >> std::ws;
        if (is.fail()) goto bad;
        else
            tmp_v.push_back((T)next_t);
      } else {
        T next_t;
        is >> next_t >> std::ws;
        if (is.fail()) goto bad;
        else
            tmp_v.push_back(next_t);
      }
    }
    is.get();  // get the final ']'.
    *v = tmp_v;  // could use std::swap to use less temporary memory, but this
    // uses less permanent memory.
  }
  if (!is.fail()) return;
 bad:
  KALDI_ERR << "ReadIntegerVector: read failure at file position "
            << is.tellg();
}

// Initialize an opened stream for writing by writing an optional binary
// header and modifying the floating-point precision.
inline void InitKaldiOutputStream(std::ostream &os, bool binary) {
  // This does not throw exceptions (does not check for errors).
  if (binary) {
    os.put('\0');
    os.put('B');
  }
  // Note, in non-binary mode we may at some point want to mess with
  // the precision a bit.
  // 7 is a bit more than the precision of float..
  if (os.precision() < 7)
    os.precision(7);
}

/// Initialize an opened stream for reading by detecting the binary header and
// setting the "binary" value appropriately.
inline bool InitKaldiInputStream(std::istream &is, bool *binary) {
  // Sets the 'binary' variable.
  // Throws exception in the very unusual situation that stream
  // starts with '\0' but not then 'B'.

  if (is.peek() == '\0') {  // seems to be binary
    is.get();
    if (is.peek() != 'B') {
      return false;
    }
    is.get();
    *binary = true;
    return true;
  } else {
    *binary = false;
    return true;
  }
}

}  // end namespace kaldi.

#endif  // KALDI_BASE_IO_FUNCS_INL_H_
