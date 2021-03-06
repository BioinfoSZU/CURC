/* Based on LzmaLib.h -- LZMA library interface
2013-01-18 : Igor Pavlov : Public domain */

#ifndef __LZMA_LIB_H
#define __LZMA_LIB_H

#include "../7zTypes.h"
#include "helper.h"
#include "VarLenDNACoder.h"
#include <vector>

using namespace std;

#define MY_STDAPI int MY_STD_CALL

#define LZMA_PROPS_SIZE 5

#ifdef DEVELOPER_BUILD
extern bool dump_after_decompression;
extern int dump_after_decompression_counter;
extern string dump_after_decompression_prefix;
#endif

const static uint8_t LZMA_CODER = 1;
const static uint8_t LZMA2_CODER = 2;
const static uint8_t PPMD7_CODER = 3;
const static uint8_t VARLEN_DNA_CODER = 11;
const static uint8_t COMPOUND_CODER_TYPE = 77;


const static int PGRC_DATAPERIODCODE_8_t = 0;
const static int PGRC_DATAPERIODCODE_16_t = 1;
const static int PGRC_DATAPERIODCODE_32_t = 2;
const static int PGRC_DATAPERIODCODE_64_t = 3;
const static int PGRC_DATAPERIODCODE_128_t = 4;

const static uint8_t PGRC_CODER_LEVEL_FAST = 1;
const static uint8_t PGRC_CODER_LEVEL_NORMAL = 2;
const static uint8_t PGRC_CODER_LEVEL_MAX = 3;

const static double COMPRESSION_ESTIMATION_UINT8_BITMAP = 0.125;
const static double COMPRESSION_ESTIMATION_BASIC_DNA = 0.250;
const static double COMPRESSION_ESTIMATION_VAR_LEN_DNA = 0.83;
const static double COMPRESSION_ESTIMATION_MIS_CNT = 0.5;
const static double COMPRESSION_ESTIMATION_MIS_SYM = 0.250;

double simpleUintCompressionEstimate(uint64_t dataMaxValue, uint64_t typeMaxValue);

char* Compress(size_t &destLen, const char *src, size_t srcLen, uint8_t coder_type, uint8_t coder_level,
        int coder_param = -1, double estimated_compression = 1);
void writeCompressed(ostream &dest, const char *src, size_t srcLen, uint8_t coder_type, uint8_t coder_level,
                     int coder_param = -1, double estimated_compression = 1);
void writeCompressed(ostream &dest, const string srcStr, uint8_t coder_type, uint8_t coder_level,
                     int coder_param = -1, double estimated_compression = 1);

char* componentCompress(ostream &dest, size_t &compLen, const char *src, size_t srcLen, uint8_t coder_type, uint8_t coder_level,
                        int coder_param = -1, double estimated_compression = 1);
void writeCompoundCompressionHeader(ostream &dest, size_t srcLen, size_t compLen, uint8_t coder_type);

void Uncompress(char* dest, size_t destLen, istream &src, size_t srcLen, uint8_t coder_type);
void Uncompress(char* dest, size_t destLen, const char* src, size_t srcLen, uint8_t coder_type);
void readCompressed(istream &src, string& dest);

template<typename T>
void readCompressed(istream &src, vector<T>& dest) {
    size_t destLen = 0;
    size_t srcLen = 0;
    uint8_t coder_type = 0;
    PgSAHelpers::readValue<uint64_t>(src, destLen, false);
    if (destLen % sizeof(T)) {
        fprintf(stderr, "Invalid output size %zu for decompressing to the vector of %zu-byte elements",
                destLen, sizeof(T));
    }
    dest.resize(destLen / sizeof(T));
    if (destLen == 0)
        return;
    PgSAHelpers::readValue<uint64_t>(src, srcLen, false);
    PgSAHelpers::readValue<uint8_t>(src, coder_type, false);
    Uncompress((char*) dest.data(), destLen, src, srcLen, coder_type);
#ifdef DEVELOPER_BUILD
    if (dump_after_decompression) {
        string dumpFileName = dump_after_decompression_prefix + (dump_after_decompression_counter < 10?"0":"");
        PgSAHelpers::writeArrayToFile(dumpFileName + PgSAHelpers::toString(dump_after_decompression_counter++),
                dest.data(), destLen);
    }
#endif
}

#endif