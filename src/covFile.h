#include "includedefine.h"

class buffer_in_chunk{
  private:
    char * compressed_buffer;
    char * buffer;

    unsigned long bufferPos;
    unsigned long bufferMax;
    uint64_t l_file_buffer;
    char * file_buffer;
    uint64_t file_bufferPos;
  public:
    buffer_in_chunk();
    ~buffer_in_chunk();
    bool BufferIsFull(unsigned int threshold = 8);
    int write(char * src, unsigned int len);
    int WriteBuffer();
    
    char * get_buffer_ptr() { return file_buffer; };
    uint64_t get_buffer_pos() { return file_bufferPos; };
};

class covReader {
  private:
    char * compressed_buffer;
    char * buffer;

    unsigned long bufferPos;
    unsigned long bufferMax;
    
    uint32_t index_begin;
    uint32_t body_begin;

    istream * IN;
    
    int IS_EOF;
    int IS_FAIL;
    size_t IS_LENGTH;
    size_t EOF_POS;
    
    std::vector<std::string> chr_names;
    std::vector<unsigned int> chr_lens;

  public:
    covReader();
    ~covReader();

    void SetInputHandle(std::istream *in_stream);
    
    // Basic functions
    int ReadBuffer();
    int read(char * dest, unsigned int len);
    int ignore(unsigned int len);
    
    bool eof();
    bool fail();

    // Input functions
    int ReadHeader();
    int FetchPos(const std::string seqname, const uint32_t start, const int strand,
      uint64_t * file_offset, uint32_t * block_start);
    int FetchRLE(const std::string seqname, const uint32_t start, const uint32_t end, const int strand,
    std::vector<int> * values, std::vector<unsigned int> * lengths);    

    // Output functions
    
    int GetChrs(std::vector<chr_entry> &chrs) {
      if(chr_names.size() > 0) {
        for(unsigned int i = 0; i < chr_names.size(); i++) {
          chrs.push_back(chr_entry(i, chr_names.at(i), chr_lens.at(i)));
        }
      }
      return(0);
    };
};

class buffer_out_chunk {
  private:
    static const int BUFFER_OUT_CAP = 65536 - 18 - 8;

    char * buffer;
    char * compressed_buffer;
    
    unsigned int buffer_pos = 0;
    unsigned int buffer_size = 0;   // number of bytes that need to be compressed
    
    unsigned int compressed_size = 0;   // number of bytes compressed that need to be written out
    
  public:
    buffer_out_chunk();
    ~buffer_out_chunk();

    unsigned int getBGZFSize() { return(compressed_size); };
    unsigned int write(char * src, unsigned int len);

    int Compress();
    int WriteToFile(ostream * OUT);
    
    unsigned int SetPos(unsigned int pos) {
      if(pos >= BUFFER_OUT_CAP) return(buffer_pos);
      buffer_pos = pos;
      if(pos > buffer_size) {
        buffer_size = pos;
      }
      return(pos);
    };
    unsigned int GetPos() { return(buffer_pos); };
    
    unsigned int write_to_pos(char * src, unsigned int len, unsigned int pos) {
      if(len + pos > BUFFER_OUT_CAP) return(0);
      SetPos(pos);
      return(write(src, len));
    };
    
    bool IsAtCap(unsigned int len) {
      if(len + buffer_pos >= BUFFER_OUT_CAP) return(true);
      return(false);
    }
};

class covWriter {
  private:
    ostream * OUT;
    
    std::vector<chr_entry> chrs;
    
    // When chrs is set, initialize these:
    std::vector< std::vector<buffer_out_chunk> > body;    // The buffers
    std::vector< std::vector<uint32_t> > block_coord_starts;    // The start coords of each bgzf
  public:
    covWriter();
    ~covWriter();
    
    void SetOutputHandle(std::ostream *out_stream);
    
    int WriteHeader(std::vector<chr_entry> chrs_to_copy);
  
    int WriteFragmentsMap(std::vector< std::pair<unsigned int, int> > * vec, 
        unsigned int chrID, unsigned int strand);
    int WriteEmptyEntry(unsigned int refID);
    int WriteEmptyEntry(unsigned int chrID, unsigned int strand) {
      return(WriteEmptyEntry(chrID + chrs.size() * strand)); };
    
    int WriteHeaderToFile();
    int WriteIndexToFile();
    int WriteToFile();
};