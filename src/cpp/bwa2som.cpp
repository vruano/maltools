#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <getopt.h>

//#define PNG_DEBUG 3
//#include <png.h>

#include "sam.h"  
#include "faidx.h"

#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include <vector>
#include <map>

using namespace std ;

#define ARROW_LENGTH 3
#define SATURATION_THRESHOLD 64

#define NUMBER_WIDTH 4
#define NUMBER_HEIGHT 6
#define NUMBER_ALIGN_HCENTER 1
#define NUMBER_ALIGN_RIGHT 2
#define NUMBER_ALIGN_VCENTER 4

#define PC_EMPTY 0
#define PC_REFERENCE 1
#define PC_SNP 2
#define PC_SINGLE 4
#define PC_INSERTION 8
#define PC_DELETION 16

typedef long postype ;

typedef struct {  
	int beg, end;  
	samfile_t *in;  
} tmpstruct_t;  

typedef vector <postype> TVI ;

char bam2char[255] ;

class TBAMreader ;

class TBAMaligner {
	public :
	TBAMaligner ( TBAMreader *_base ) ;
	void set_range () ;
	bool draw_single_read ( const bam1_t *b, void *data ) ;
	bool draw_paired_read ( const bam1_t *b, void *data) ;
	postype get_start () { return start ; }

	postype single_offset ;
	string chr ;
	bool allow_all_snps , allow_all_indels ;
	int indel_variance ;

	protected :
	bool paint_single_read ( const bam1_t *b, void *data , int y ) ;
	bool paint_single_read_cigar ( const bam1_t *b, void *data , int y ) ;
	
	TBAMreader *base ;
	postype start , end , size ;
} ;

class TBAMreader {
	public :
	TBAMreader () {} ;
	void init ( string _bamfile , string _region , int _mapq = 0 ) ;
	void set_options ( string options ) ;
	void set_region ( string s ) { region = s ; }
	
	// BAM methods
	static int fetch_func(const bam1_t *b, void *data) ;
	static int fetch_func2(const bam1_t *b, void *data) ;
	void read_bam_file ( int mode ) ;

	// PNG methods
	inline postype get_start () { return tmp.beg ; }
	inline postype get_end () { return tmp.end ; }

	char *refseq ;
	string refseq_file , bam_file_name ;
	bamFile bam_out_file ;
	bool fix_single_reads ;
	
	TBAMaligner *draw ;
	bool o_single , o_pairs , o_arrows , o_snps , o_faceaway , o_inversions , o_linkpairs , o_colordepth ;
	bool o_noscale , o_readqual , o_text ;
	int total_snps , total_reads ;
	int highlight_from , highlight_to ;
	bool use_highlight ;
	
	private :

	void abort_(const char * s, ...) ;

	// BAM variables
	string bam_file , region ;
	tmpstruct_t tmp;
	int mapq ;
	
} ;

TBAMreader bam_aligner ; // Neccessary for hack around BAM needing static function

char IUPAC[256] ;
map < string , map < int , char > > snps ;
map < string , map < int , map < string , bool > > > insertions , deletions ;
map < string , int > read_counter ;


//////////////////////////////////////////////////////////////////////////////////////
// TBAMaligner
TBAMaligner::TBAMaligner ( TBAMreader *_base ) {
	base = _base ;
	base->draw = this ;
	single_offset = 0 ;
	allow_all_snps = false ;
	allow_all_indels = false ;
	indel_variance = 0 ;
}


void TBAMaligner::set_range () {
	start = base->get_start() + 1 ;
	end = base->get_end() + 1 ;
	size = end - start ;
}


bool TBAMaligner::draw_single_read ( const bam1_t *b, void *data ) {
	return paint_single_read ( b , data , 0 ) ; //b->core.pos % single_offset ) ;
}

bool TBAMaligner::draw_paired_read ( const bam1_t *b, void *data) {
	// Paint read
	return paint_single_read ( b , data , abs(b->core.isize) ) ;
}

bool TBAMaligner::paint_single_read ( const bam1_t *b, void *data , int y ) {
	return paint_single_read_cigar ( b , data , y ) ; // QUICK HACK!
/*	
	if ( y < vstart || y >= vend ) return ;

	if ( b->core.n_cigar > 1 ) {
		paint_single_read_cigar ( bucket , b , data , y ) ;
		return ;
	}
	
	postype from = b->core.pos ;
	postype to = from + b->core.l_qseq - 1 ;
	

	
	// SNPs
	if ( base->o_snps && base->refseq ) {
		uint8_t *s = bam1_seq ( b ) ;
		postype p = b->core.pos + 1 ;
		for ( postype a = 0 ; a < b->core.l_qseq ; a++ , p++ ) {
			if ( p < start ) continue ;
			if ( p >= end ) break ;

			if ( *(base->refseq+p-start) != bam2char[bam1_seqi(s,a)] ) {
			}
		}
	}

	*/
}

bool TBAMaligner::paint_single_read_cigar ( const bam1_t *b, void *data , int y ) {
	bool ret = true ;
	uint32_t *cigdata = bam1_cigar(b) ;
	postype p = b->core.pos ;
	int rp = 0 ;
	uint8_t *s = bam1_seq ( b ) ;
	
	for ( int cigcnt = 0 ; ret && cigcnt < b->core.n_cigar ; cigcnt++ ) {
		uint32_t ciglen = cigdata[cigcnt] >> 4 ;
		uint32_t cigtype = cigdata[cigcnt] & 15 ;
		
		if ( cigtype == BAM_CMATCH ) {
			for ( int b = 0 ; b < ciglen ; b++ , p++ , rp++ ) {
				char rc = bam2char[bam1_seqi(s,rp)] ;
				if ( *(base->refseq+p-start+1) != rc ) {
					int gpos = p-start+2 ; // Why is this "+2" and not "+1"? It works, but why? WHY??? FIXME
					if ( rc != 'N' && rc != 'n' && snps[chr].find(gpos) != snps[chr].end() && ( IUPAC[rc] & snps[chr][gpos] ) ) {
//						cout << gpos << ":" << (char)rc << " : OK" << endl ;
					} else {
//						cout << gpos << ":" << (char)rc << " : BAD" << endl ;
						ret = false ;
					}
					
				} else if ( rc == 'N' || rc == 'n' ) {
					ret = false ;
				} else {
				}
			}
			if ( allow_all_snps ) ret = true ;
		} else if ( cigtype == BAM_CINS ) {
			string seq ;
			int gpos = p-start+1 ; // Here it's "+1" - see above
			for ( int a = 0 ; a < ciglen ; a++ ) seq += bam2char[bam1_seqi(s,rp+a)] ;
//			cout << chr << ":" << gpos << "\t" << seq << endl ;
			if ( allow_all_indels ) ret = true ;
			else {
				bool again = true ;
				for ( int v = -indel_variance ; again && ret && v <= indel_variance ; v++ ) {
					if ( insertions.find(chr) != insertions.end() && insertions[chr].find(gpos) != insertions[chr].end() && insertions[chr][gpos].find(seq) != insertions[chr][gpos].end() ) again = false ;
//					else ret = false ;
				}
				ret = !again ;
			}
			rp += ciglen ;

		} else if ( cigtype == BAM_CDEL ) {
			string seq ;
			int gpos = p-start+1 ; // Here it's "+1" - see above
			for ( int b = 0 ; b < ciglen ; b++ , p++ ) {
				seq += *(base->refseq+p-start+1) ;
			}
//			cout << chr << ":" << gpos << "\t" << seq << endl ;
			if ( allow_all_indels ) ret = true ;
			else {
				bool again = true ;
				for ( int v = -indel_variance ; again && ret && v <= indel_variance ; v++ ) {
					if ( deletions.find(chr) != deletions.end() && deletions[chr].find(gpos) != deletions[chr].end() && deletions[chr][gpos].find(seq) != deletions[chr][gpos].end() ) again = false ;
//					else ret = false ;
				}
				ret = !again ;
			}

		} else if ( cigtype == BAM_CREF_SKIP ) { // UNTESTED
			ret = false ;
			rp += ciglen ;
		} else if ( cigtype == BAM_CSOFT_CLIP ) { // UNTESTED
			ret = false ;
			rp += ciglen ;
		} else if ( cigtype == BAM_CHARD_CLIP ) { // UNTESTED
			ret = false ;
			p += ciglen ;
		} else if ( cigtype == BAM_CPAD ) { // UNTESTED
			ret = false ;
			p += ciglen ;
		} else ret = false ;
		
	}
	
	return ret ;
}






//////////////////////////////////////////////////////////////////////////////////////
// TBAMreader

void TBAMreader::init ( string _bam_file , string _region , int _mapq ) 
{
	fix_single_reads = false ;
	use_highlight = false ;
	bam_file = _bam_file ;
//	sam_out_file = NULL ;
	bam_out_file = NULL ;
	region = _region ;
	mapq = _mapq ;
	draw = NULL ;
	refseq = NULL ;
	o_single = o_pairs = o_arrows = o_snps = o_faceaway = o_inversions = o_linkpairs = o_colordepth = false ;
	o_noscale = o_readqual = o_text = false ;
}

void TBAMreader::set_options ( string options ) {
	const char *last = options.c_str() ;
	char *c = (char*) options.c_str() ;
	vector <string> ov ;
	for ( c++ ; *c ; c++ ) {
		if ( *c == ',' ) {
			*c = 0 ;
			ov.push_back ( last ) ;
			last = c + 1 ;
		}
	}
	ov.push_back ( last ) ;
	for ( int a = 0 ; a < ov.size() ; a++ ) {
		if ( ov[a] == "snps" ) o_snps = true ;
		else if ( ov[a] == "pairs" ) o_pairs = true ;
		else if ( ov[a] == "arrows" ) o_arrows = true ;
		else if ( ov[a] == "single" ) o_single = true ;
		else if ( ov[a] == "faceaway" ) o_faceaway = true ;
		else if ( ov[a] == "inversions" ) o_inversions = true ;
		else if ( ov[a] == "linkpairs" ) o_linkpairs = true ;
		else if ( ov[a] == "colordepth" ) o_colordepth = true ;
		else if ( ov[a] == "noscale" ) o_noscale = true ;
		else if ( ov[a] == "readqual" ) o_readqual = true ;
		else if ( ov[a] == "text" ) o_text = true ;
	}
	
	if ( o_single ) draw->single_offset = 50 ;
}


void TBAMreader::read_bam_file ( int mode ) {
	if ( mode == 0 && draw == NULL ) abort_ ( "No analysis class instanced!\n" ) ;

	
	tmp.beg = 0 ;
	tmp.end = 0;   
	tmp.in = samopen(bam_file.c_str(), "rb", 0); 


	// WTF does this loop do???
//	cout << tmp.in->header->text << endl ;
	for ( char *c = tmp.in->header->text ; *c ; c++ ) {
		if ( *c == '\n' && *(c+1) == '@' && *(c+2) == 'C' && *(c+3) == 'O' ) {
			string s ;
			use_highlight = false ;
			for ( char *d = c+5 ; *d > 13 ; d++ ) {
				if ( *d == ' ' ) {
					if ( s == "HIGHLIGHT" ) use_highlight = true ;
					s = "" ;
				} else s += *d ;
			}
			if ( use_highlight ) {
				highlight_from = atoi ( (char*) s.c_str() ) ;
				const char *d ;
				for ( d = s.c_str() ; *(d-1) != '-' ; d++ ) ;
				highlight_to = atoi ( d ) ;
			}
		}
	}
	
	
	int ref;  
	bam_index_t *idx;  
	bam_plbuf_t *buf;  
	idx = bam_index_load(bam_file.c_str()); // load BAM index  
	if (idx == 0) {  
		fprintf(stderr, "BAM indexing file is not available.\n");  
		return ;  
	}  
	bam_parse_region(tmp.in->header, region.c_str(), &ref,  &tmp.beg, &tmp.end); // parse the region  
	
	if ( bam_aligner.draw ) bam_aligner.draw->chr = region ;
	
	if ( !refseq_file.empty() ) {
		faidx_t *fai = fai_load ( refseq_file.c_str() ) ;
		int len = tmp.end - tmp.beg + 2 ;
		refseq = fai_fetch ( fai , region.c_str() , &len ) ;
		if ( len < tmp.end - tmp.beg ) tmp.end = len + tmp.beg ;
/*		if ( tmp.beg == 0 ) { // No range set
			tmp.beg = 1 ;
			tmp.end = len ;
		}*/
		fai_destroy ( fai ) ;
	}
	
//	if ( width * 50 < tmp.end - tmp.beg ) o_arrows = false ;

        cerr << "opening file....";
        cerr <<  bam_file_name;	
        cerr << "\n";
	if ( bam_out_file == NULL ) {
		bam_out_file = bam_open ( bam_file_name.c_str() , "wh" ) ;
                cerr << "done\n";
		bam_header_write ( bam_out_file , tmp.in->header ) ;
	}
        cerr << "continue\n";
	
	if ( bam_aligner.draw ) draw->set_range () ;
	total_snps = 0 ;
	total_reads = 0 ;
	if ( mode == 0 ) bam_fetch(tmp.in->x.bam, idx, ref, tmp.beg, tmp.end, NULL, fetch_func);  
	else if ( mode == 1 ) bam_fetch(tmp.in->x.bam, idx, ref, tmp.beg, tmp.end, NULL, fetch_func2);  
	bam_index_destroy(idx);  

//	cout << "Total SNPs : " << total_snps << endl ;
}

int TBAMreader::fetch_func2(const bam1_t *b, void *data) {
	map <string,int>::iterator i = read_counter.find ( string ( bam1_qname(b) ) ) ;
	if ( i != read_counter.end() ) {
		bam1_t *b2 = (bam1_t*) b ; // AAAARGH MY EYES!!!
		if ( b2->core.flag & BAM_FPROPER_PAIR ) b2->core.flag ^= BAM_FPROPER_PAIR ;
		b2->core.flag |= BAM_FMUNMAP ;
		read_counter.erase ( i ) ;
	}
	
	bam_write1 ( bam_aligner.bam_out_file , b ) ;
	return 0 ;
}

int TBAMreader::fetch_func(const bam1_t *b, void *data) {  
	if ( b->core.qual < bam_aligner.mapq ) return 0 ;

	if ( b->core.qual <= 0 ) return 0 ; // Minimum alignment quality

	bam_aligner.total_reads++ ;
	bool ok = false ;
	if ( b->core.flag & BAM_FPROPER_PAIR ) {
		ok = bam_aligner.draw->draw_paired_read ( b , data ) ;
	} else if ( b->core.flag & BAM_FUNMAP ) {
	} else if ( b->core.flag & BAM_FMUNMAP ) {
		ok = bam_aligner.draw->draw_single_read ( b , data ) ;
	} else if ( b->core.isize != 0 ) {
		ok = bam_aligner.draw->draw_paired_read ( b , data ) ;
	} else {
		ok = bam_aligner.draw->draw_single_read ( b , data ) ;
	}
	
	if ( ok ) {
		
		if ( bam_aligner.fix_single_reads && ( ((b->core.flag&BAM_FPROPER_PAIR)>0) || !((b->core.flag&BAM_FMUNMAP)>0) ) ) {
			string qname ( bam1_qname(b) ) ;
			read_counter[qname]++ ;
			if ( read_counter[qname] == 2 ) read_counter.erase ( qname ) ;
		}
		
//		samwrite ( bam_aligner.sam_out_file , b ) ;
		bam_write1 ( bam_aligner.bam_out_file , b ) ;
	}
	
	return 0;  
}  



void TBAMreader::abort_(const char * s, ...)
{
	va_list args;
	va_start(args, s);
	vfprintf(stderr, s, args);
	fprintf(stderr, "\n");
	va_end(args);
	abort();
}



// MAIN STUFF

void init_iupac () {
	for ( int a = 0 ; a < 256 ; a++ ) IUPAC[a] = 0 ;
	IUPAC['A'] = 1 ;
	IUPAC['C'] = 2 ;
	IUPAC['G'] = 4 ;
	IUPAC['T'] = 8 ;
	IUPAC['R'] = IUPAC['A'] | IUPAC['G'] ;
	IUPAC['Y'] = IUPAC['C'] | IUPAC['T'] ;
	IUPAC['S'] = IUPAC['G'] | IUPAC['C'] ;
	IUPAC['W'] = IUPAC['A'] | IUPAC['T'] ;
	IUPAC['K'] = IUPAC['G'] | IUPAC['T'] ;
	IUPAC['M'] = IUPAC['A'] | IUPAC['C'] ;
	IUPAC['B'] = IUPAC['C'] | IUPAC['G'] | IUPAC['T'] ;
	IUPAC['D'] = IUPAC['A'] | IUPAC['G'] | IUPAC['T'] ;
	IUPAC['H'] = IUPAC['A'] | IUPAC['C'] | IUPAC['T'] ;
	IUPAC['V'] = IUPAC['A'] | IUPAC['C'] | IUPAC['G'] ;
	IUPAC['N'] = IUPAC['A'] | IUPAC['C'] | IUPAC['G'] | IUPAC['T'] ;
	
	for ( int a = 'A' ; a <= 'Z' ; a++ ) IUPAC[a-'A'+'a'] = IUPAC[a] ;
}

void read_indel_list ( string file ) {
	cout << "Reading " << file << " ... " ;
	
	ifstream in ( file.c_str() ) ;
	while ( !in.eof() ) {
		string s ;
		getline ( in , s ) ;

		string chr ;
		int pos ;
		char v ;
		string seq ;
		const char *c = s.c_str() ;
		if ( *c == 0 ) break ; // This is the end
		while ( *c != 9 ) chr += *c++ ;
		c++ ;
		pos = atoi ( c ) ;
		while ( *c != 9 ) c++ ;
		c++ ;
		v = *c++ ;
		while ( *c && *c != 9 ) seq += *c++ ;
		
		if ( chr == "Chr" ) continue ;
		
		if ( v == '+' ) insertions[chr][pos][seq] = true ;
		else deletions[chr][pos][seq] = true ;
		
//		cout << chr << ":" << pos << " " << (char)v << " : " << seq << endl ;
	}

	cout << "done.\n" ;
	in.close() ;
}

void read_snp_list ( string file ) {
	cout << "Reading " << file << " ... " ;
	
	ifstream in ( file.c_str() ) ;
	while ( !in.eof() ) {
		string s ;
		getline ( in , s ) ;

		string chr ;
		char ref , alt ;
		int pos ;
		const char *c = s.c_str() ;
		if ( *c == 0 ) break ; // This is the end
		while ( *c != 9 ) chr += *c++ ;
		c++ ;
		pos = atoi ( c ) ;
		while ( *c != 9 ) c++ ;
		c++ ;
		ref = *c ;
		c += 2 ;
		alt = *c ;
		
		snps[chr][pos] = IUPAC[ref] | IUPAC[alt] ;
		
//		cout << chr << ":" << pos << " " << (char) ref << "/" << (char) alt << endl ;
	}
	
	cout << "done.\n" ;
	in.close() ;
}

/*
// OBSOLETE
void do_fix_single_reads ( string infile , string outfile ) {
	cout << "Fixing single reads... " << outfile << endl ;
	ifstream in ( infile.c_str() ) ;
	ofstream out ( outfile.c_str() ) ;
	
	while ( !in.eof() ) {
		string s , readname ;
		getline ( in , s ) ;
		if ( s.empty() ) continue ;
		
		if ( s[0] == '@' ) {
			out << s << endl ;
			continue ;
		}
		
		const char *c ;
		for ( c = s.c_str() ; *c && *c != 9 ; readname += *c++ ) ;
		if ( !*c ) { out << s << endl ; continue ; }
		
		if ( read_counter.find ( readname ) == read_counter.end() ) {
			out << s << endl ;
			continue ;
		}
		
		c++ ;
		int key = atoi ( c ) ;
		while ( *c && *c != 9 ) c++ ;
		if ( !*c ) { out << s << endl ; continue ; }
		if ( key & BAM_FPROPER_PAIR ) key ^= BAM_FPROPER_PAIR ;
		key |= BAM_FMUNMAP ;
		out << readname << "\t" << key << c << endl ;
		read_counter.erase ( readname ) ;
	}
	
	in.close() ;
	out.close() ;
	remove ( outfile.c_str() ) ;
}
*/

int die_usage () {
	cout << "Usage : bwa2som --bam=FILE --ref=FILE --out=OUT_BAM_FILE <options>" << endl ;
	cout << "-S --snps    FILE     SNP list [optional]" << endl ;
	cout << "-I --indels  FILE     InDel list [optional]" << endl ;
	cout << "-R --region  REGION   Region [optional]" << endl ;
	cout << "-f --fix              Fix single-read flag (doubles runtime, for LookSeq etc.) [optional]" << endl ;
	cout << "-i --iv      INT      InDel position variance [optional]" << endl ;
	return 1 ;
}

void run_regions ( string region , string ref_file , int mode ) {
	if ( region.empty() ) { // TODO check for existence of .fai file!
		string file = ref_file + ".fai" ;
		ifstream in ( file.c_str() ) ;
		while ( !in.eof() ) {
			string s , chr ;
			getline ( in , s ) ;
			for ( const char *c = s.c_str() ; *c && *c != 9 && *c != 32 ; chr += *c++ ) ;
			if ( chr.empty() ) continue ;
			bam_aligner.set_region ( chr ) ;
			bam_aligner.read_bam_file ( mode ) ;
		}
	} else {
		bam_aligner.read_bam_file ( mode ) ;
	}
	bam_close ( bam_aligner.bam_out_file ) ;
	bam_index_build ( bam_aligner.bam_file_name.c_str() ) ;
}

int main(int argc, char **argv) {
	for ( int a = 0 ; a < 256 ; a++ ) bam2char[a] = '?' ;
	bam2char[1] = 'A' ;
	bam2char[2] = 'C' ;
	bam2char[4] = 'G' ;
	bam2char[8] = 'T' ;
	bam2char[15] = 'N' ;
	
	init_iupac () ;

	string view = "indel" ;
	string bam_file , ref_file , region , bam_file_name ;
        string tmp_dir = "";
	string snp_list , indel_list ;
	int indel_variance = 0 ;
	int mapq = 0 ;
	bool fixbam = false ;
	static struct option long_options[] = {
                { "tmpdir" , optional_argument, 0, 't' },
		{ "bam" , optional_argument , 0 , 'b' } ,
		{ "ref" , optional_argument , 0 , 'r' } ,
		{ "region" , optional_argument , 0 , 'R' } ,
		{ "mapq" , optional_argument , 0 , 'm' } ,
		{ "snps" , optional_argument , 0 , 'S' } ,
		{ "indels" , optional_argument , 0 , 'I' } ,
		{ "out" , optional_argument , 0 , 'o' } ,
		{ "fix" , optional_argument , 0 , 'f' } ,
		{ "iv" , optional_argument , 0 , 'i' } ,
		{ 0 , 0 , 0 , 0 }
	} ;
	int c = 0 ;
	while ( -1 != ( c = getopt_long (argc, argv, "",long_options, NULL) ) ) {
		switch ( c ) {
			case 0 : break ;
			case 'b' : bam_file = optarg ; break ;
			case 'r' : ref_file = optarg ; break ;
			case 'R' : region = optarg ; break ;
			case 'S' : snp_list = optarg ; break ;
			case 'I' : indel_list = optarg ; break ;
			case 'o' : bam_file_name = optarg ; break ;
			case 'f' : fixbam = true ; break ;
			case 'm' : mapq = atoi ( optarg ) ; break ;
                        case 't' : tmp_dir = optarg; break;
			case 'i' : indel_variance = atoi ( optarg ) ; break ;
		}
	}
	
	if ( bam_file.empty() ) return die_usage () ;
	if ( bam_file_name.empty() ) return die_usage () ;

	if ( !snp_list.empty() && snp_list != "all" ) read_snp_list ( snp_list ) ;
	if ( !indel_list.empty() && indel_list != "all" ) read_indel_list ( indel_list ) ;

	bam_aligner.refseq_file = ref_file ;
	bam_aligner.init( bam_file , region , mapq ) ;

	bam_aligner.fix_single_reads = fixbam ;

	if ( bam_aligner.fix_single_reads ) {
		if ( 1 ) { // Use /tmp directory
			string templt = (tmp_dir.empty()) ? "bwa2somXXXXXX" : tmp_dir + "/bwa2somXXXXXX";
  			char *t = new char[templt.size() + 1];
			strcpy(t,templt.c_str());
			int fd;
			fd = mkstemp(t);
			close ( fd ) ;
			bam_aligner.bam_file_name = t ;
                        delete[] t;
		} else {
			bam_aligner.bam_file_name = bam_file_name + ".tmp" ;
		}

		cerr << "Using temp file " << bam_aligner.bam_file_name << endl ;
	} else bam_aligner.bam_file_name = bam_file_name ;


	TBAMaligner aligner ( &bam_aligner ) ;
	aligner.indel_variance = indel_variance ;
	if ( snp_list == "all" ) aligner.allow_all_snps = true ;
	if ( indel_list == "all" ) aligner.allow_all_indels = true ;

	run_regions ( region , ref_file , 0 ) ;
	
	if ( bam_aligner.fix_single_reads ) {
		cerr << "Fixing BAM file ..." << endl ;
		bam_file = bam_aligner.bam_file_name ;
		
		bam_aligner = TBAMreader () ;
		bam_aligner.refseq_file = ref_file ;
		bam_aligner.init( bam_file , region , mapq ) ;
		bam_aligner.bam_file_name = bam_file_name ;
		run_regions ( region , ref_file , 1 ) ;

		remove ( bam_file.c_str() ) ;
		bam_file += ".bai" ;
		remove ( bam_file.c_str() ) ;
	}

	return 0;
}

/*
\rm bwa2som test.sam* test.bam ; g++ bwa2som.cpp -O3 -o bwa2som -lz -L . -lbam ; time ./bwa2som --ref=/nfs/team112/refseq/plasmodium/falciparum/3D7_pm.fa --bam=/nfs/users/nfs_m/mm6/ftp/pf.bwa/bam/PP0002-C.bam --snps=/nfs/team112/snplists/plasmodium/falciparum/bwa.976K.2010-07-26.snps --indels=/nfs/team112/pipeline/development/alignment_and_snp_calling/bwa.pf/deli/indel.data.list --sam=test.sam ; samtools view -S -u test.sam | samtools sort -n -o - - | samtools fixmate - - | samtools sort - test


\rm bwa2som test.* ; icpc bwa2som.cpp -O3 -o bwa2som -lz -L . -lbam ; time ./bwa2som --ref=/nfs/team112/refseq/plasmodium/falciparum/3D7_pm.fa --bam=/nfs/users/nfs_m/mm6/ftp/pf.bwa/bam/PP0002-C.bam --snps=/nfs/team112/snplists/plasmodium/falciparum/bwa.976K.2010-07-26.snps --indels=/nfs/team112/pipeline/development/alignment_and_snp_calling/bwa.pf/deli/indel.data.list --out=test.bam ; samtools pileup -vcf /nfs/team112/refseq/plasmodium/falciparum/3D7_pm.fa test.bam > test.pileup ; gawk '{ if ( $3 == "*" ) print $0 }' test.pileup | wc -l


\rm bwa2som test.* ; icpc bwa2som.cpp -O3 -o bwa2som -lz -L . -lbam ; time ./bwa2som --ref=/nfs/team112/refseq/plasmodium/falciparum/3D7_pm.fa --bam=/nfs/users/nfs_m/mm6/ftp/pf.bwa/bam/PH0003-C.bam --snps=/nfs/team112/snplists/plasmodium/falciparum/bwa.976K.2010-07-26.snps --out=PH0003-C.filtered.noindels.bam

samtools pileup -cf /nfs/team112/refseq/plasmodium/falciparum/3D7_pm.fa test.bam | ./bam2snps --snps=/nfs/team112/snplists/plasmodium/falciparum/bwa.976K.2010-07-26.snps --minq=27 > test.snps

samtools mpileup -Augf /nfs/team112/refseq/plasmodium/falciparum/3D7_pm.fa test.bam | bcftools view -bvcg - | bcftools view -





\rm bwa2som ; g++ bwa2som.cpp -O3 -o bwa2som -lpng -L . -lbam ; time ./bwa2som --bam=/nfs/users/nfs_m/mm6/ftp/ag/bam/AC0001-C.bam --options=snps,pairs,arrows,single,faceaway,inversions,linkpairs,colordepth --ref=/nfs/users/nfs_m/mm6/ftp/ag/Anopheles_gambiae.clean.fa --region="2L:1-200000" --png=2L.a.png


\rm bwa2som ; g++ bwa2som.cpp -O3 -o bwa2som -lpng -L . -lbam ; time ./bwa2som --bam="ftp://ftp.sanger.ac.uk/pub/team112/ag/bam/AC0001-C.bam" --options=pairs,arrows,single,faceaway,inversions,colordepth,snps --ref=/nfs/users/nfs_m/mm6/ftp/ag/Anopheles_gambiae.clean.fa --region="2L" --png=2L.a.png

\rm bwa2som ; g++ bwa2som.cpp -O3 -o bwa2som -lpng -L . -lbam ; cp bwa2som ~/wwwdev_data_marker3/..

*/
