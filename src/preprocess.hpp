/**
Software License:
-----------------

Copyright (c) 2021 Zexuan Zhu<zhuzx@szu.edu.cn>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

PgRC:
-----------------

The PgRC is an in-memory algorithm for compressing the DNA stream of FASTQ
datasets, based on the idea of building an approximation of the shortest
common superstring over high-quality reads. CURC is based on the architecture
of PgRC and also uses parts of PgRC codes in backend encoding
(mainly variable-length encoding and pairing data encoding).

The PgRC copyright is as follows:

Copyright (c) 2020 Tomasz M. Kowalski, Szymon Grabowski All Rights Reserved.

See also the PgRC web site:
  https://github.com/kowallus/PgRC for more information.

*/

#ifndef CURC_PREPROCESS_HPP
#define CURC_PREPROCESS_HPP

#include "Param.hpp"
void compress(Param& );

#endif //CURC_PREPROCESS_HPP
