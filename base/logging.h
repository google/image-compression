// Copyright 2015 Google Inc. All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef BASE_LOGGING_H_
#define BASE_LOGGING_H_

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>

class LogNull {
 public:
  LogNull() {}
  ~LogNull() {}
  std::ostream& GetStream() {
    static std::ostream null_stream(NULL);
    return null_stream;
  }
};

#ifdef _DEBUG

class LogFatal {
 public:
  LogFatal(const char* file, int line) { str_ << file << ":" << line << ": "; }
  ~LogFatal() {
    const std::string s = str_.str();
    std::cerr << s << std::endl;
    abort();
  }
  std::ostream& GetStream() { return str_; }

 private:
  std::ostringstream str_;
};

#define LOG_IF_NOT(condition)                             \
  !(condition) ? LogFatal(__FILE__, __LINE__).GetStream() \
               : LogNull().GetStream()

#define DCHECK_OP(a, b, op) LOG_IF_NOT((a)op(b))

#define DCHECK_LE(a, b) DCHECK_OP((a), (b), <=)
#define DCHECK_GE(a, b) DCHECK_OP((a), (b), >=)
#define DCHECK_LT(a, b) DCHECK_OP((a), (b), <)
#define DCHECK_GT(a, b) DCHECK_OP((a), (b), >)
#define DCHECK_EQ(a, b) DCHECK_OP((a), (b), ==)
#define DCHECK_NE(a, b) DCHECK_OP((a), (b), !=)
#define DCHECK(a) LOG_IF_NOT((a))

#else

#define DCHECK(condition) LogNull().GetStream()
#define DCHECK_EQ(a, b) LogNull().GetStream()
#define DCHECK_NE(a, b) LogNull().GetStream()
#define DCHECK_GT(a, b) LogNull().GetStream()
#define DCHECK_LT(a, b) LogNull().GetStream()
#define DCHECK_GE(a, b) LogNull().GetStream()
#define DCHECK_LE(a, b) LogNull().GetStream()

#endif

#endif  // BASE_LOGGING_H_
