#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <getopt.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>

using namespace std ;

#define BASE_A 1
#define BASE_C 2
#define BASE_G 4
#define BASE_T 8

typedef map < string , map < int , map < string , bool > > > vmap ;

char IUPAC[256] ;
map < string , map < int , char > > snps ;
vmap insertions , deletions ;
map < string , bool > chromosomes ;
FILE *snpout = NULL ;
FILE *indelout = NULL ;

void init_iupac () {
	for ( int a = 0 ; a < 256 ; a++ ) IUPAC[a] = 0 ;
	IUPAC['A'] = BASE_A ;
	IUPAC['C'] = BASE_C ;
	IUPAC['G'] = BASE_G ;
	IUPAC['T'] = BASE_T ;
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

void read_snp_list ( string file ) {
//	cout << "Reading " << file << " ... " ;
	
	int cnt = 0 ;
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
		if ( chromosomes.find(chr) == chromosomes.end() ) chromosomes[chr] = true ;
		cnt++ ;
		
//		cout << chr << ":" << pos << " " << (char) ref << "/" << (char) alt << endl ;
	}
	
//	cout << cnt << " done.\n" ;
	in.close() ;
}

void read_indel_list ( string file ) {
//	cout << "Reading InDel file " << file << " ... " ;
	
	int cnt = 0 ;
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

		if ( chromosomes.find(chr) == chromosomes.end() ) chromosomes[chr] = true ;
		cnt++ ;
		
//		cout << chr << ":" << pos << " " << (char)v << " : " << seq << endl ;
	}

//	cout << cnt << " done.\n" ;
	in.close() ;
}

bool is_listed_indel ( string chr , int pos , string var ) {
	if ( var == "*" ) return true ;
	string seq = var.substr ( 1 ) ;
	vmap *v ;
	if ( var[0] == '+' ) {
		v = &insertions ;
	} else if ( var[0] == '-' ) {
		v = &deletions ;
	} else return false ;
	if ( (*v).find(chr) == (*v).end() ) return false ;
	if ( (*v)[chr].find(pos) == (*v)[chr].end() ) return false ;
	if ( (*v)[chr][pos].find(seq) == (*v)[chr][pos].end() ) return false ;
	return true ;
}

void parse_mpileup ( int minq ) { // FIXME!!!!!! Works for SNPs only!
	while ( !cin.eof() ) {
		string s ;
		getline ( cin , s ) ;
		string chr , ref , bases , quals ;
		int pos , cov ;
		const char *c = s.c_str() ;
		if ( !*c ) continue ; // Blank line
		while ( *c && *c != 9 ) chr += *c++ ;
		if ( *c ) c++ ; else continue ;
		pos = atoi ( c ) ;
		
		if ( chromosomes.find ( chr ) == chromosomes.end() ) continue ; // No chromosome we care about!
		
		
		while ( *c && *c != 9 ) c++ ;
		if ( *c ) c++ ; else continue ;
		while ( *c && *c != 9 ) ref += *c++ ;
		if ( *c ) c++ ; else continue ;
		cov = atoi ( c ) ;
		while ( *c && *c != 9 ) c++ ;
		if ( *c ) c++ ; else continue ;
		while ( *c && *c != 9 ) bases += *c++ ;
		if ( *c ) c++ ; else continue ;
		while ( *c && *c != 9 ) quals += *c++ ;

		if ( ref.length() != 1 ) continue ; // Paranoia
		
		if ( snpout && ref != "*" && snps[chr].find ( pos ) != snps[chr].end() ) {
			char refb = ref[0] ;
			
			int cnt[256] ;
			cnt['A'] = cnt['C'] = cnt['G'] = cnt['T'] = 0 ;
			
			const char *b = bases.c_str() ;
			const char *q = quals.c_str() ;
			while ( *b ) {
				char rb , rq ;
				rb = rq = '?' ;
				if ( *b == '^' || *b == '$' ) { // Begin/end of read, skip
					if ( *b == '$' ) {
						b++ ;
					} else {
						b += 2 ;
					}
					continue ;
				}
				if ( *b == '.' || *b == ',' ) { // Reference
					if ( !*q ) break ; // Paranoia
					rb = refb ;
					rq = *q++ ;
					b++ ;
				} else if ( *b >= 'A' && *b <= 'Z' ) {
					if ( !*q ) break ; // Paranoia
					rb = *b++ ;
					rq = *q++ ;
				} else if ( *b >= 'a' && *b <= 'z' ) {
					if ( !*q ) break ; // Paranoia
					rb = (*b++) - 'a' + 'A' ;
					rq = *q++ ;
				} else if ( *b == '+' || *b == '-' ) {
					char type = *b++ ;
					int len = atoi ( b ) ;
					while ( *b >= '0' && *b <= '9' ) b++ ;
					string seq ;
					for ( int a = 0 ; a < len ; a++ ) seq += *b++ ;
					// TODO store indel
					continue ;
				} else if ( *b == '*' ) { // Deletion, skip
					if ( !*q ) break ; // Paranoia
					b++ ;
					q++ ;
					continue ;
				} else {
					cout << "ERROR" << endl << s << endl << bases << endl << "!!!! : " << b << endl ;
					exit ( 0 ) ;
				}
				
				if ( rb == '?' ) {
					cerr << "ERROR\n" << s << endl ;
					exit ( 0 ) ;
				}
				
				rq -= 33 ;
				if ( rq < minq ) continue ; // Quality filter
				
				cnt[rb]++ ;
				
			}
			
			if ( *b || *q ) {
				cout << "ERROR" << endl << s << endl ;
				if ( *b ) cout << "SEQ : " << b << endl ;
				if ( *q ) cout << "QUA : " << q << endl ;
				exit ( 0 ) ;
			}
			
			if ( ( ( BASE_A & IUPAC[snps[chr][pos]] ) && cnt['A'] > 0 ) ||
				 ( ( BASE_C & IUPAC[snps[chr][pos]] ) && cnt['C'] > 0 ) ||
				 ( ( BASE_G & IUPAC[snps[chr][pos]] ) && cnt['G'] > 0 ) ||
				 ( ( BASE_T & IUPAC[snps[chr][pos]] ) && cnt['T'] > 0 ) ) {
				cerr << "ERROR : Bad base" << endl << s << endl ;
				exit ( 0 ) ;
			}
	
			if ( cnt['A'] > 0 ) fprintf ( snpout , "%7d\t%s\t%d\t%c\t%c\n" , cnt['A'] , chr.c_str() , pos , 'A' , refb ) ;
			if ( cnt['C'] > 0 ) fprintf ( snpout , "%7d\t%s\t%d\t%c\t%c\n" , cnt['C'] , chr.c_str() , pos , 'C' , refb ) ;
			if ( cnt['G'] > 0 ) fprintf ( snpout , "%7d\t%s\t%d\t%c\t%c\n" , cnt['G'] , chr.c_str() , pos , 'G' , refb ) ;
			if ( cnt['T'] > 0 ) fprintf ( snpout , "%7d\t%s\t%d\t%c\t%c\n" , cnt['T'] , chr.c_str() , pos , 'T' , refb ) ;
//		} else if ( indelout && ref == "*" && ( insertions[chr].find ( pos ) != insertions[chr].end() || deletions[chr].find ( pos ) != deletions[chr].end() ) ) {
//			fprintf ( indelout , "%s\n" , s.c_str() ) ;
		}
	}
}


void parse_pileup ( int minq ) {
	int snpcount = 0 ;
	int indelcount = 0 ;
	while ( !cin.eof() ) {
		string s ;
		getline ( cin , s ) ;
		string chr , ref , bases , quals ;
		int pos , cov ;
		const char *c = s.c_str() ;
		if ( !*c ) continue ; // Blank line
		while ( *c && *c != 9 ) chr += *c++ ;
		if ( *c ) c++ ; else continue ;
		pos = atoi ( c ) ;
		
		if ( chromosomes.find ( chr ) == chromosomes.end() ) continue ; // No chromosome we care about!
		
		
		while ( *c && *c != 9 ) c++ ;
		if ( *c ) c++ ; else continue ;
		while ( *c && *c != 9 ) ref += *c++ ;
		if ( *c ) c++ ; else continue ;
		
		for ( int u = 0 ; u < 4 ; u++ ) {
			while ( *c && *c != 9 ) *c++ ;
			c++ ;
		}
		
		cov = atoi ( c ) ;
		while ( *c && *c != 9 ) c++ ;
		if ( *c ) c++ ; else continue ;

		if ( ref.length() != 1 ) continue ; // Paranoia
		
		if ( snpout && ref != "*" && snps[chr].find ( pos ) != snps[chr].end() ) { // SNP position

			while ( *c && *c != 9 ) bases += *c++ ;
			if ( *c ) c++ ; else continue ;
			while ( *c && *c != 9 ) quals += *c++ ;

			char refb = ref[0] ;
			
			int cnt[256] ;
			cnt['A'] = cnt['C'] = cnt['G'] = cnt['T'] = 0 ;
			snpcount++ ;
			
			const char *b = bases.c_str() ;
			const char *q = quals.c_str() ;
			while ( *b ) {
				char rb , rq ;
				rb = rq = '?' ;
				if ( *b == '^' || *b == '$' ) { // Begin/end of read, skip
					if ( *b == '$' ) {
						b++ ;
					} else {
						b += 2 ;
					}
					continue ;
				}
				if ( *b == '.' || *b == ',' ) { // Reference
					if ( !*q ) break ; // Paranoia
					rb = refb ;
					rq = *q++ ;
					b++ ;
				} else if ( *b >= 'A' && *b <= 'Z' ) {
					if ( !*q ) break ; // Paranoia
					rb = *b++ ;
					rq = *q++ ;
				} else if ( *b >= 'a' && *b <= 'z' ) {
					if ( !*q ) break ; // Paranoia
					rb = (*b++) - 'a' + 'A' ;
					rq = *q++ ;
				} else if ( *b == '+' || *b == '-' ) {
					char type = *b++ ;
					int len = atoi ( b ) ;
					while ( *b >= '0' && *b <= '9' ) b++ ;
					string seq ;
					for ( int a = 0 ; a < len ; a++ ) seq += *b++ ;
					// TODO store indel
					continue ;
				} else if ( *b == '*' ) { // Deletion, skip
					if ( !*q ) break ; // Paranoia
					b++ ;
					q++ ;
					continue ;
				} else {
					cout << "ERROR" << endl << s << endl << bases << endl << "!!!! : " << b << endl ;
					exit ( 0 ) ;
				}
				
				if ( rb == '?' ) {
					cerr << "ERROR\n" << s << endl ;
					exit ( 0 ) ;
				}
				
				rq -= 33 ;
				if ( rq < minq ) continue ; // Quality filter
				
				cnt[rb]++ ;
				
			}
			
			if ( *b || *q ) {
				cout << "ERROR" << endl << s << endl ;
				if ( *b ) cout << "SEQ : " << b << endl ;
				if ( *q ) cout << "QUA : " << q << endl ;
				exit ( 0 ) ;
			}
			
			if ( ( ( BASE_A & IUPAC[snps[chr][pos]] ) && cnt['A'] > 0 ) ||
				 ( ( BASE_C & IUPAC[snps[chr][pos]] ) && cnt['C'] > 0 ) ||
				 ( ( BASE_G & IUPAC[snps[chr][pos]] ) && cnt['G'] > 0 ) ||
				 ( ( BASE_T & IUPAC[snps[chr][pos]] ) && cnt['T'] > 0 ) ) {
				cerr << "ERROR : Bad base" << endl << s << endl ;
				exit ( 0 ) ;
			}
	
			if ( cnt['A'] > 0 ) fprintf ( snpout , "%7d\t%s\t%d\t%c\t%c\n" , cnt['A'] , chr.c_str() , pos , 'A' , refb ) ;
			if ( cnt['C'] > 0 ) fprintf ( snpout , "%7d\t%s\t%d\t%c\t%c\n" , cnt['C'] , chr.c_str() , pos , 'C' , refb ) ;
			if ( cnt['G'] > 0 ) fprintf ( snpout , "%7d\t%s\t%d\t%c\t%c\n" , cnt['G'] , chr.c_str() , pos , 'G' , refb ) ;
			if ( cnt['T'] > 0 ) fprintf ( snpout , "%7d\t%s\t%d\t%c\t%c\n" , cnt['T'] , chr.c_str() , pos , 'T' , refb ) ;
			
		} else if ( indelout && ref == "*" && ( insertions[chr].find ( pos ) != insertions[chr].end() || deletions[chr].find ( pos ) != deletions[chr].end() ) ) { // InDel call
			string v1 , v2 ;
			int cnt1 , cnt2 ;
			
			while ( *c && *c != 9 ) v1 += *c++ ;
			if ( *c ) c++ ; else continue ;
			while ( *c && *c != 9 ) v2 += *c++ ;
			if ( *c ) c++ ; else continue ;

			cnt1 = atoi ( c ) ;
			while ( *c && *c != 9 ) c++ ;
			if ( *c ) c++ ; else continue ;
			cnt2 = atoi ( c ) ;
			while ( *c && *c != 9 ) c++ ;
			if ( *c ) c++ ; else continue ;

			if ( is_listed_indel ( chr , pos , v1 ) && cnt1 > 0 ) fprintf ( indelout , "%7d\t%s\t%d\t%s\t%s\n" , cnt1 , chr.c_str() , pos , v1.c_str() , "*" ) ;
			if ( is_listed_indel ( chr , pos , v2 ) && cnt2 > 0 ) fprintf ( indelout , "%7d\t%s\t%d\t%s\t%s\n" , cnt2 , chr.c_str() , pos , v2.c_str() , "*" ) ;
			indelcount++ ;
		} else if ( indelout && ref != "*" && ( insertions[chr].find ( pos ) != insertions[chr].end() || deletions[chr].find ( pos ) != deletions[chr].end() ) ) { // InDel position with no call
			fprintf ( indelout , "%7d\t%s\t%d\t%s\t%s\n" , cov , chr.c_str() , pos , "*" , "*" ) ;
			indelcount++ ;
		}
	}
//	cout << "snpcount\t" << snpcount << endl ;
//	cout << "indelcount\t" << indelcount << endl ;
}


int main(int argc, char **argv) {
	string snpfile , indelfile ;
	string mode = "pileup" ;
	int minq = 0 ;
	static struct option long_options[] = {
		{ "mode" , optional_argument , 0 , 'm' } ,
		{ "minq" , optional_argument , 0 , 'q' } ,
		{ "snplist" , optional_argument , 0 , 's' } ,
		{ "indellist" , optional_argument , 0 , 'i' } ,
		{ "snpout" , optional_argument , 0 , 'S' } ,
		{ "indelout" , optional_argument , 0 , 'I' } ,
		{ 0 , 0 , 0 , 0 }
	} ;
	int c = 0 ;
	while ( -1 != ( c = getopt_long (argc, argv, "",long_options, NULL) ) ) {
		switch ( c ) {
			case 0 : break ;
			case 'm' : mode = optarg ; break ;
			case 's' : snpfile = optarg ; break ;
			case 'i' : indelfile = optarg ; break ;
			case 'S' : snpout = fopen ( optarg , "w" ) ; break ;
			case 'I' : indelout = fopen ( optarg , "w" ) ; break ;
			case 'q' : minq = atoi(optarg) ; break ;
		}
	}
	
	if ( snpfile.empty() ) {
		cerr << "USAGE :" << endl << "samtools pileup -cf REFSEQ BAMFILE | bam2snps [--snplist=SNPFILE] [--indellist=INDELFILE] [--snpout=SNP_OUT_FILE] [--indelout=INDEL_OUT_FILE] [--minq=MINIMUM_SNP_BASE_QUALITY]" << endl ;
		cerr << "OR    :" << endl << "samtools mpileup -Bf REFSEQ BAMFILE | bam2snps --mode=mpileup [--snplist=SNPFILE] [--snpout=SNP_OUT_FILE] [--minq=MINIMUM_SNP_BASE_QUALITY]" << endl ;
		exit ( 0 ) ;
	}
	
	init_iupac () ;
	if ( !snpfile.empty() ) read_snp_list ( snpfile ) ;
	if ( !indelfile.empty() ) read_indel_list ( indelfile ) ;

	if ( mode == "pileup" ) parse_pileup ( minq ) ;
	else if ( mode == "mpileup" ) parse_mpileup ( minq ) ;
	else { cerr << "Unknown mode " << mode << endl ; return 1 ; }
	
	if ( snpout ) fclose ( snpout ) ;
	if ( indelout ) fclose ( indelout ) ;
	
	return 0 ;
}


/*
\rm bam2snps ; icpc -O3 bam2snps.cpp -o bam2snps ; samtools pileup -vcf /nfs/team112/refseq/plasmodium/falciparum/3D7_pm.fa test.bam | ./bam2snps --snplist=/nfs/team112/snplists/plasmodium/falciparum/bwa.976K.2010-07-26.snps --indellist=/nfs/team112/pipeline/development/alignment_and_snp_calling/bwa.pf/deli/indel.data.list --minq=27 --snpout=test.snps --indelout=test.indels

\rm bam2snps ; icpc -O3 bam2snps.cpp -o bam2snps
samtools pileup -Bf /nfs/team112/refseq/plasmodium/falciparum/3D7_pm.fa BAMFILE | ./bam2snps
*/
