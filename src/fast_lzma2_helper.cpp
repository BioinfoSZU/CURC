#include <bits/stdc++.h>
#include "fast_lzma2_helper.hpp"
#include "../third_party/fast-lzma2/fast-lzma2.h"

namespace lzma2 {

    static void exit_fail(const char *msg)
    {
        fputs(msg, stderr);
        exit(1);
    }

    static int compress_file(FL2_CStream *fcs, FILE *fin, FILE *fout)
    {
        unsigned char in_buffer[8 * 1024];
        unsigned char out_buffer[4 * 1024];
        FL2_inBuffer in_buf = { in_buffer, sizeof(in_buffer), sizeof(in_buffer) };
        FL2_outBuffer out_buf = { out_buffer, sizeof(out_buffer), 0 };
        size_t res = 0;
        size_t in_size = 0;
        size_t out_size = 0;
        do {
            if (in_buf.pos == in_buf.size) {
                in_buf.size = fread(in_buffer, 1, sizeof(in_buffer), fin);
                in_size += in_buf.size;
                in_buf.pos = 0;
            }
            res = FL2_compressStream(fcs, &out_buf, &in_buf);
            if (FL2_isError(res))
                goto error_out;

            fwrite(out_buf.dst, 1, out_buf.pos, fout);
            out_size += out_buf.pos;
            out_buf.pos = 0;

        } while (in_buf.size == sizeof(in_buffer));
        do {
            res = FL2_endStream(fcs, &out_buf);
            if (FL2_isError(res))
                goto error_out;

            fwrite(out_buf.dst, 1, out_buf.pos, fout);
            out_size += out_buf.pos;
            out_buf.pos = 0;
        } while (res);

        return 0;

        error_out:
        fprintf(stderr, "Error: %s\n", FL2_getErrorName(res));
        return 1;
    }

    static int decompress_file(FL2_DStream *fds, FILE *fin, FILE *fout)
    {
        unsigned char in_buffer[4 * 1024];
        unsigned char out_buffer[8 * 1024];
        FL2_inBuffer in_buf = { in_buffer, sizeof(in_buffer), sizeof(in_buffer) };
        FL2_outBuffer out_buf = { out_buffer, sizeof(out_buffer), 0 };
        size_t res;
        size_t in_size = 0;
        size_t out_size = 0;
        do {
            if (in_buf.pos == in_buf.size) {
                in_buf.size = fread(in_buffer, 1, sizeof(in_buffer), fin);
                in_size += in_buf.size;
                in_buf.pos = 0;
            }
            res = FL2_decompressStream(fds, &out_buf, &in_buf);
            if (FL2_isError(res))
                goto error_out;

            fwrite(out_buf.dst, 1, out_buf.pos, fout);
            out_size += out_buf.pos;
            out_buf.pos = 0;
        } while (res && in_buf.size);

        return 0;

        error_out:
        fprintf(stderr, "Error: %s\n", FL2_getErrorName(res));
        return 1;
    }

    static void open_files(const char *infile, const char *outfile, FILE **fin, FILE **fout)
    {
        *fin = fopen(infile, "rb");
        if (*fin == NULL)
            exit_fail("Cannot open input file.\n");

        *fout = fopen(outfile, "wb");
        if (*fout == NULL)
            exit_fail("Cannot open output file.\n");
    }

    static void create_init_fl2_stream_compression(FL2_CStream **fcs, int level, int thread_num)
    {
        *fcs = FL2_createCStreamMt(thread_num, 0);
        if (*fcs == nullptr)
            exit_fail("Cannot allocate compression context.\n");

        size_t res = FL2_initCStream(*fcs, level);
        if (FL2_isError(res)) {
            fprintf(stderr, "Error: %s\n", FL2_getErrorName(res));
            exit(1);
        }
    }

    static void create_init_fl2_stream_decompression(FL2_DStream **fds)
    {
        *fds = FL2_createDStreamMt(16);
        if (*fds == nullptr)
            exit_fail("Cannot allocate decompression context.\n");

        size_t res = FL2_initDStream(*fds);
        if (FL2_isError(res)) {
            fprintf(stderr, "Error: %s\n", FL2_getErrorName(res));
            exit(1);
        }
    }

    void lzma2_compress(const char *infile, const char *outfile, int level, int thread_num) {
        FILE *fin;
        FILE *fout;
        FL2_CStream *fcs;
        create_init_fl2_stream_compression(&fcs, level, thread_num);
        open_files(infile, outfile, &fin, &fout);
        compress_file(fcs, fin, fout);
        fclose(fout);
        fclose(fin);
        FL2_freeCStream(fcs);
    }

    void lzma2_decompress(const char *infile, const char *outfile) {
        FILE *fin;
        FILE *fout;
        FL2_DStream *fds;
        create_init_fl2_stream_decompression(&fds);
        open_files(infile, outfile, &fin, &fout);
        decompress_file(fds, fin, fout);
        fclose(fout);
        fclose(fin);
        FL2_freeDStream(fds);
    }

} // namespace lzma2

