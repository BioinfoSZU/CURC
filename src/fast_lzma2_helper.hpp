#ifndef CURC_FAST_LZMA_HELPER_HPP
#define CURC_FAST_LZMA_HELPER_HPP

namespace lzma2 {

    void lzma2_compress(const char *infile, const char *outfile, int level, int thread_num);

    void lzma2_decompress(const char *infile, const char *outfile);

} // namespace lzma2

#endif //CURC_FAST_LZMA_HELPER_HPP
