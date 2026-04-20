#if !defined(BROTLI_BASE64_LUT_)
#define BROTLI_BASE64_LUT_
static const uint8_t kIsBase64[256] = {
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,  /* 43 '+', 47 '/' (45 '-' is 0) */
    1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,  /* 48-57 '0'-'9' (61 '=' is 0) */
    0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,  /* 65-79 'A'-'O' */
    1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,  /* 80-90 'P'-'Z' */
    0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,  /* 97-111 'a'-'o' */
    1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,  /* 112-122 'p'-'z' */
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
};
#endif

#define HistogramType FN(Histogram)

static void FN(BuildAndStoreEntropyCodes)(MemoryManager* m, BlockEncoder* self,
    const HistogramType* histograms, const size_t histograms_size,
    const size_t alphabet_size, HuffmanTree* tree,
    const uint8_t* is_base64_histogram,
    size_t* storage_ix, uint8_t* storage) {
  const size_t table_size = histograms_size * self->histogram_length_;
  self->depths_ = BROTLI_ALLOC(m, uint8_t, table_size);
  self->bits_ = BROTLI_ALLOC(m, uint16_t, table_size);
  if (BROTLI_IS_OOM(m)) return;

  {
    size_t i;
    for (i = 0; i < histograms_size; ++i) {
      size_t ix = i * self->histogram_length_;
      if (self->histogram_length_ == 256 && is_base64_histogram && is_base64_histogram[i]) {
        size_t k;
        memset(&self->depths_[ix], 0, 256);
        for (k = 0; k < 256; ++k) {
          if (kIsBase64[k]) {
            self->depths_[ix + k] = 6;
          }
        }
        BrotliConvertBitDepthsToSymbols(&self->depths_[ix], 256, &self->bits_[ix]);
        BrotliStoreHuffmanTree(&self->depths_[ix], 256, tree, storage_ix, storage);
      } else {
        BuildAndStoreHuffmanTree(&histograms[i].data_[0], self->histogram_length_,
            alphabet_size, tree, &self->depths_[ix], &self->bits_[ix],
            storage_ix, storage);
      }
    }
  }
}

#undef HistogramType
