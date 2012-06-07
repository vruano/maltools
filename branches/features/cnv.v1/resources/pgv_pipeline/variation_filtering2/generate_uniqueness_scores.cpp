#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <iostream>

using namespace std ;

typedef string index_type ;
typedef uint32_t count_type ;
typedef map <string,string> ssmap ;
typedef map <index_type,count_type> simap ;
typedef map <string,vector<uint8_t> > mvimap ;

/// Loads an entire FASTA file into a vector of TChromosome.
void load_genome ( string filename , ssmap &genome ) {
    fstream filestr;
    filestr.open ( filename.c_str() , fstream::in);
    string name ;
    while ( !filestr.eof() ) {
        string s ;
        getline ( filestr , s ) ;
        if ( s[0] == '>' ) {
            const char *c = s.c_str() ;
            while ( *c == '>' || *c == ' ' ) c++ ;
            name = c ;
        } else genome[name] += s ;
    }
    filestr.close() ;
    
    for ( ssmap::iterator i = genome.begin() ; i != genome.end() ; i++ ) {
        for ( uint32_t p = 0 ; p < i->second.length() ; p++ ) {
            switch ( i->second[p] ) {
                case 'A' :
                case 'C' :
                case 'G' :
                case 'T' : break ;
                case 'a' : i->second[p] = 'A' ; break ;
                case 'c' : i->second[p] = 'C' ; break ;
                case 'g' : i->second[p] = 'G' ; break ;
                case 't' : i->second[p] = 'T' ; break ;
                default : i->second[p] = 'N' ; break ;
            }
        }
    }
}

index_type get_index ( const char *start , int len ) {
    string s1 ( start , len ) ;
    string s2 ( len , 'N' ) ;
    for ( int p = 0 ; p < len ; p++ ) {
        switch ( s1[p] ) {
            case 'A' : s2[len-p-1] = 'T' ; break ;
            case 'C' : s2[len-p-1] = 'G' ; break ;
            case 'G' : s2[len-p-1] = 'C' ; break ;
            case 'T' : s2[len-p-1] = 'A' ; break ;
        }
    }
    return s1 > s2 ? s1 : s2 ;
}

int main ( int argc , char *argv[] ) {
    int index_start = 10 ;
    int index_end = 50 ;
    int index_step = 1 ;
//    string reference = "sources_test/reference.fasta" ;
    string reference = argv[1] ;
    ssmap genome ;
    load_genome ( reference , genome ) ;

    mvimap min ;
    for ( ssmap::iterator i = genome.begin() ; i != genome.end() ; i++ ) {
        min[i->first].resize(i->second.length(),0) ;
    }

    for ( int index_size = index_start ; index_size <= index_end ; index_size += index_step ) {
        cerr << "INDEX " << index_size << endl ;
        
        // Create index
        simap index ;
        for ( ssmap::iterator i = genome.begin() ; i != genome.end() ; i++ ) {
//            cout << i->first << endl ;
            uint32_t chrlen = i->second.length() ;
            const char *base = i->second.c_str() ;
            for ( uint32_t pos = 0 ; pos < chrlen - index_size ; pos++ ) {
                index[get_index(base+pos,index_size)]++ ;
            }
        }
        
        // Check index
        for ( ssmap::iterator i = genome.begin() ; i != genome.end() ; i++ ) {
            int32_t chrlen = i->second.length() ;
            const char *base = i->second.c_str() ;
            for ( int32_t pos = 0 ; pos < chrlen ; pos++ ) {
                if ( min[i->first][pos] > 0 ) continue ; // Had that
                
                bool pure = true ;
                for ( int offset = -index_size ; offset <= 0 ; offset++ ) {
                    if ( pos + offset < 0 || pos + offset + index_size >= chrlen ) continue ;
                    if ( index[get_index(base+pos+offset,index_size)] == 1 ) continue ;
                    pure = false ;
                    break ;
                }
                
                if ( pure ) min[i->first][pos] = index_size ;
            }
        }
        
    }

    // Output
    for ( ssmap::iterator i = genome.begin() ; i != genome.end() ; i++ ) {
        int32_t chrlen = i->second.length() ;
        for ( int32_t pos = 0 ; pos < chrlen ; pos++ ) {
            int v = min[i->first][pos] ;
            if ( v == 0 ) v = 99 ;
            cout << i->first << "\t" << ( pos+1 ) << "\t" << v << endl ;
        }
    }

    return 0 ;
}


// \rm generate_uniqueness_scores ; icpc generate_uniqueness_scores.cpp -o generate_uniqueness_scores ; ./generate_uniqueness_scores sources_test/reference.fasta > uniqueness.scores
// 10-50 at step 1 runs ~4h
