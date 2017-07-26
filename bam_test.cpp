#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>

#include "htslib/sam.h"

using std::cout;
using std::cerr;
using std::cin;
using std::endl;
using std::string;
using std::ofstream;
using std::to_string;


/*
KHASH_MAP_INIT_INT(bin, bins_t)                                                 
typedef khash_t(bin) bidx_t;  

typedef struct {                                                                
  int32_t n, m;                                                               
  uint64_t *offset;                                                           
} lidx_t; 

struct __hts_idx_t {                                                            
  int fmt, min_shift, n_lvls, n_bins;                                         
  uint32_t l_meta;                                                            
  int32_t n, m;                                                               
  uint64_t n_no_coor;                                                         
  bidx_t **bidx;                                                              
  lidx_t *lidx;                                                               
  uint8_t *meta; // MUST have a terminating NUL on the end                    
  struct {                                                                    
    uint32_t last_bin, save_bin;                                            
    int last_coor, last_tid, save_tid, finished;                            
    uint64_t last_off, save_off;                                            
    uint64_t off_beg, off_end;                                              
    uint64_t n_mapped, n_unmapped;                                          
  } z; // keep internal states                                                
};                                            

typedef struct {
    uint32_t read_rest:1, finished:1, is_cram:1, dummy:29;
    int tid, beg, end, n_off, i;
    int curr_tid, curr_beg, curr_end;
    uint64_t curr_off;
    hts_pair64_t *off;
    hts_readrec_func *readrec;
    struct {
        int n, m;
        int *a;
    } bins;
} hts_itr_t;
*/
int main(int argc, char **argv) {

  htsFile *hts_fp;
  if ((hts_fp = sam_open(argv[1], "r")) == NULL) {
    cerr << "Can't open sam/bam file:" << argv[1] << endl;
    exit(EXIT_FAILURE);
  }
  hts_idx_t *hts_idx = sam_index_load(hts_fp, argv[2]);
  if (hts_idx == NULL) {
    cerr << "Can't get index for bam file" << endl;
    exit(EXIT_FAILURE);
  }
  
  bam_hdr_t * bam_hdr = bam_hdr_read(hts_fp->fp.bgzf);

  // hts_idx_t *hts_idx = bam_index_load(argv[2]);
  int n_job = 10;
  // int total_job = bam_hdr->target_len[0] - 1;
  cout << bam_hdr->target_len[0] << endl;
  int total_job = 100000;
  int avg_job = total_job / n_job;
#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < n_job; ++i) { 
    htsFile *hts_fp_inner;
    if ((hts_fp_inner = sam_open(argv[1], "r")) == NULL) {
      cerr << "Can't open sam/bam file:" << argv[1] << endl;
      exit(EXIT_FAILURE);
    }
    

    cout << i << endl;
    int tid = 0;
    int beg = 9500000 + avg_job * i, end;
    //int beg = avg_job * i, end;
    if ((end = beg + avg_job) > (4800000 + 4800000)) {
      end = total_job;
    }
    cout << beg << ' ' << end << endl;
    // hts_itr_t *hts_itr = bam_itr_queryi(hts_idx, tid, beg, end);
    hts_itr_t *hts_itr = sam_itr_queryi(hts_idx, tid, beg, end);
    bam1_t *bam1 = bam_init1();
    ofstream ofs;
    string fn = argv[3] + to_string(i);
    ofs.open(fn.c_str(), ofstream::out);
    int n_read = 0;
    while ((n_read = sam_itr_next(hts_fp_inner, hts_itr, bam1)) >= 0) {
#if 0
      ofs << bam1->l_data << ' ' << bam1->m_data << endl;
      for (uint32_t j = 0; j < bam1->l_data; ++j) {
        ofs << bam1->data[j] << endl;  
      }
      ofs << bam1->id << endl;
      ofs << bam1->core.tid << endl;
      ofs << bam1->core.pos << endl;
      ofs << bam1->core.bin << endl;
      ofs << bam1->core.qual << endl;
      ofs << bam1->core.l_qname << endl;
      ofs << bam1->core.flag << endl;
      ofs << bam1->core.unused1 << endl;
      ofs << bam1->core.l_extranul << endl;
      ofs << bam1->core.n_cigar << endl;
      ofs << bam1->core.l_qseq << endl;
      ofs << bam1->core.mtid << endl;
      ofs << bam1->core.mpos << endl;
      ofs << bam1->core.isize << endl;
#endif
      uint8_t *q = bam_get_seq(bam1);
      string qseq;
      for(int i=0; i< bam1->core.l_qseq ; i++){
              qseq.push_back(seq_nt16_str[bam_seqi(q,i)]); //gets nucleotide id and converts them into IUPAC id.
      }
      ofs << qseq << endl;
    }
    ofs.close();
    sam_itr_destroy(hts_itr);
    bam_destroy1(bam1);
    sam_close(hts_fp_inner);
  }
  bam_hdr_destroy(bam_hdr);
  sam_close(hts_fp);
  return 0;
}
/*

typedef struct {
    int32_t tid;
    int32_t pos;
    uint16_t bin;
    uint8_t qual;
    uint8_t l_qname;
    uint16_t flag;
    uint8_t unused1;
    uint8_t l_extranul;
    uint32_t n_cigar;
    int32_t l_qseq;
    int32_t mtid;
    int32_t mpos;
    int32_t isize;
} bam1_core_t;

typedef struct {
    bam1_core_t core;
    int l_data;
    uint32_t m_data;
    uint8_t *data;
#ifndef BAM_NO_ID
    uint64_t id;
#endif
} bam1_t;
*/
