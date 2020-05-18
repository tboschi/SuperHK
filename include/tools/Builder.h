#ifndef _Builder
#define _Builder

#include <vector>
#include <map>
#include <cstdlib>
#include <string>
#include <iostream>

#include "TypeWrapper.h"

#include "TTree.h"
#include "TArrayF.h"
#include "TString.h"


class AlignFriend;

class memory_map{

  public:   

  memory_map()
  {
     ndims = 0;
     dims[0] = 0;
     dims[1] = 0;
     dims[2] = 0;
     dims[3] = 0;
     dims[4] = 0;
     address = 0;
  }

  void * address;
  int  size;
  std::string title;
  std::string Type;

  // how many dimensions entries per dimension
  int  dims[5];  
  int  ndims;

};


class Builder 
{

 public:

  Builder( TTree * );
 ~Builder();

  void  Get( const char *, void * ptr );
  void* Get( const char * ); 
  void  GetEntry( unsigned int x );
  int   GetEntries() { return _ttree->GetEntries(); }
  unsigned GetCurrentEntry() { return fCurrentEntry; }
  

  void SetBranch( const char * );
  void SetAllBranches(int verbose=0);
  void EnableUsedListOnly( );

  std::vector<TString> GetListOfKeys();
  std::vector<TString> GetListOfBranchTitles();
  int   GetSize( const char * );

  bool Existence( TString key ) { return ( _branch_space.count( key.Data() ) > 0 ? true :false ); }

  void Get( TString , TypeWrapper<int>&      );
  void Get( TString , TypeWrapper<unsigned>& );
  void Get( TString , TypeWrapper<float>&    );
  void Get( TString , TypeWrapper<double>&   );

  TArrayF * GetTArrayF( TString );

  void AddToUseList( TString key ) { use_list[key] = 1 ;}

  void AlignmentTest( AlignFriend * , unsigned skip = 10 );

  void SetVerbosity( int v ) { kVerbosity = v ; }

protected: 
  TTree * _ttree;
  TList * ListOfFriends;
  std::vector< int > LongFriendList;

  void SetAllBranches(TTree * , int verbose=0);
  std::map< TString , int > use_list;
  std::map< TString , TArrayF * >                    _tarrayf_data;

  std::map< std::string , memory_map* >              _branch_space;
  std::map< std::string , int         >              _var_size;

  std::map< TString , TypeWrapper <int> *       >    _int_data;
  std::map< TString , TypeWrapper <float> *     >    _float_data;
  std::map< TString , TypeWrapper <double> *    >    _double_data;
  std::map< TString , TypeWrapper <unsigned> *  >    _unsigned_data;
  std::map< TString , TypeWrapper <char> *      >    _void_data;
  std::map< TString , TypeWrapper <ULong64_t> * >    _ulong64_data;
  std::map< TString , TypeWrapper <Long64_t> *  >    _long64_data;


  unsigned fCurrentEntry;
  bool kBranchesSet ;
  
  int kVerbosity ;
};





#endif



