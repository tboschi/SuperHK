#include <string>
#include <sstream>

#include "Builder.h"
#include "TLeaf.h"
#include "TBranch.h"
#include "TObjArray.h"
#include "TList.h"
#include "TFriendElement.h"
#include "TChain.h"
#include "TFile.h"

#include "AlignFriend.h"

#include "toolsLibVer.h"

Builder::Builder( TTree * tree )
{

  kVerbosity = 0 ;
  kBranchesSet = false ;
  _ttree = tree;
  ListOfFriends = NULL;
  ListOfFriends = _ttree->GetListOfFriends();
 
  // if _ttree is a TChain, the perform a first check to 
  // see if the friends have more 
  // entries than the number of entries of the first tree 
  // in the chain. This is done so that the correct entry 
  // in the friend can be called when the user specifies the 
  // global entry in the chain. 
  //
  // If _ttree is a TTree then all of the friends should
  // at least as many entries as _ttree but there is no ambiguity
  // arising from a global entry vs. a local (to a chain element) 
  // entry. 
  //
  // See Builder::GetEntry
  // 
  if( ListOfFriends && _ttree->IsA() == TChain::Class() )
  if(0)
  {

      TTree * ltree;
      TFriendElement * fe;
      _ttree->LoadTree(0);
      for( int i = 0 ; i < ListOfFriends->GetSize() ; i++ )
      {
         fe = (TFriendElement*) ListOfFriends->At(i);
         ltree = fe->GetTree();
         //if( ltree->GetEntries() > _ttree->GetEntries() )
         if( ltree->GetEntries() > _ttree->GetTree()->GetEntries() )
         {
             std::cout << " FriendEntries: " << ltree->GetEntries() << " Parent: " << _ttree->GetTree()->GetEntries() << std::endl;
             LongFriendList.push_back(1);
         }   
         else 
             LongFriendList.push_back(0);
      }

  }// end of check for friends

  _ttree->GetEntry(0);
 
 toolsLibVer * tlvp = toolsLibVer::GetLibrary();
 LibraryVerMaster::GetMaster()->AddLibrary( tlvp ); 
}

Builder::~Builder(  )
{
   std::map< std::string , memory_map*>::iterator _i;
   for( _i = _branch_space.begin() ; _i != _branch_space.end() ; _i++ )
   {
      free( _i->second->address );
      delete _i->second;
   }

}


void* Builder::Get( const char * branch  )
{
   std::string branchname( branch );

   int size;
   if( _branch_space.count( branchname ) ) 		//similar to index checking
      return _branch_space[ branchname ]->address ;
   else
      std::cout << " Warning in Builder::Get - branch -" << branch << "- was not found. " << std::endl;

   return 0;
}


void Builder::Get( const char * branch, void * ptr  )
{
   std::string branchname( branch );
   
   int size;
   if( _branch_space.count( branchname ) ) 
   {
      std::cout <<  " --- " << branchname << " " << ptr << " " << _branch_space[ branchname ]->address << std::endl;
      size = _branch_space[ branchname ]->size;
      memcpy( ptr, _branch_space[ branchname ]->address, size );

   }
   else
      std::cout << " Warning in Builder::Get - branch " << branch << " was not found. " << std::endl;

}

void Builder::SetBranch( const char * name   )	//using map to save branches, seems faster while prcessing data?
{
   TBranch * branch = 0;
   std::string ltitle;
   std::string * _str;
   std::string   sub;
   std::string   Type;
   std::string   branchname( name );
   std::string::size_type loc0;
   std::string::size_type loc1;
   int size  = 0 ;  // defaults to non-array variables

   int type_size;
   int total_size;

   struct memory_map * _store = new struct memory_map; 

   branch = _ttree->GetBranch( name );
   if( ! branch ) { std::cerr << "Builder::SetBranch Branch: " << branchname << " does not exist. " << std::endl; abort();} 

   ltitle =  branch->GetTitle() ;  

   _str = new std::string( branch->GetTitle()  );
   TString classname = branch->GetClassName();

   if ( classname == "TArrayF" )
   {
       _tarrayf_data[ branchname ] = new TArrayF; 
       //std::cout<< " --- Storing TArrayF " << branchname << " " << std::endl;  

      _store->address = &_tarrayf_data[ branchname ];
      _store->size    = sizeof( TArrayF* );
      _store->title   =  ltitle ;
      _store->Type    = "T";
      _store->dims[0] = 1;
      _store->ndims   = 1;

      _branch_space[ branchname ] = _store;
      _ttree->SetBranchAddress( branchname.c_str(),  &_tarrayf_data[ branchname ] );
      return;
   }


   total_size = 1; // the default is to _not_ be an array
   loc0  = _str->find("[");
   loc1  = _str->find("]");
   int count = 0;
   if( loc0 == std::string::npos )
   {
      // single variables are just 1-dimensional arrays with a single entry
      _store->dims[0] = 1;
      _store->ndims   = 1;
   }
   
   if ( kVerbosity >= 2)
      std::cout << "Builder::SetBranch Loading string: " << _str->c_str() << " " << _store->ndims << " " << std::endl;
   while( loc0 != std::string::npos && loc1 != std::string::npos )
   {
      sub   = _str->substr( loc0+1 , loc1-loc0-1 ); 
      size  = atoi( sub.c_str() );

      if ( kVerbosity >= 2)
        std::cout << "Builder::SetBranch sub: " << sub << " " << size <<  " "  << count <<  " " << _store << std::endl;
         
      // for variable length array indices controlled by another leaf ( ie nring )    
      if( size < 1 ) 
      {
       if( ! _var_size.count(sub) )
       {
         TLeaf * ll = _ttree->GetLeaf( sub.c_str() );

         // this is an error condition for the case myVar[0] or myVar[] 
         if( ! ll ) size = 1 ;
         else       size = (int) ll->GetMaximum();

         if ( size <= 1 ) size = (int) _ttree->GetMaximum( sub.c_str() );
         if ( size == 0 ) size = 1 ; 

         _var_size[ sub ] = size ;
         if ( kVerbosity >= 2)
           std::cout << "Builder::SetBranch Size of variable : " << sub << " " << size << " at finish " <<  branchname.c_str() << std::endl; 
       }
       else 
       {	
         size = _var_size[ sub ] ; 
         if ( kVerbosity >= 2)
            std::cout << "Builder::SetBranch Using size " << size << " for var " << sub << " from var_size map " << std::endl;
       }	
      }
      
      if ( kVerbosity >= 2)
        std::cout << "Builder::SetBranch sub: " << sub << " size: " << size << " tot: " << total_size << std::endl;

      total_size *= size;

      loc0  = _str->find("[",loc0+1);
      loc1  = _str->find("]",loc0+1);

      _store->ndims = count + 1 ; 
      _store->dims[count++] = size;
   }

   // now search for the data type:
   loc0 = _str->find("/");
   Type = _str->substr(loc0+1);

   type_size = 0; 
   if( Type == "D" ) type_size = sizeof( double );
   if( Type == "F" ) type_size = sizeof( float );
   if( Type == "I" ) type_size = sizeof( int );
   if( Type == "i" ) type_size = sizeof( unsigned int );
   if( Type == "l" ) type_size = sizeof( ULong64_t );
   if( Type == "L" ) type_size = sizeof( Long64_t  );
   if( Type == "C" ) type_size = 256 * sizeof( char );
   if( type_size == 0 ) return;  // unsupported type

   total_size *= type_size;
   if( total_size == 0 ) { std::cerr << "Builder::SetBranch Size of branch : " << branchname << " computed to be zero. abort. " << std::endl; abort();} 

   _store->address = calloc( total_size , sizeof(char) );         

   if ( kVerbosity >= 2)
     std::cout << "Builder::SetBranch Address: " << _store->address  << std::endl;

   _store->size    = total_size;
   _store->title   = ltitle ;
   _store->Type    = Type;  
  

   if( Type == "D" ) 
      _double_data[ branchname ] = 
            new TypeWrapper<double>( (double*)_store->address, _store->dims, _store->ndims, branchname  );

   if( Type == "F" ) 
      _float_data[ branchname ]  = 
            new TypeWrapper<Float_t>( (float*)_store->address, _store->dims, _store->ndims , branchname );

   if( Type == "I" ) 
      _int_data[ branchname ]    =  
            new TypeWrapper<Int_t>( (int*)_store->address, _store->dims, _store->ndims , branchname );

   if( Type == "i" ) 
      _unsigned_data[ branchname ] = 
            new TypeWrapper<unsigned>( (unsigned*)_store->address, _store->dims, _store->ndims, branchname  );

   if( Type == "C" ) 
      _void_data[ branchname ]  = 
            new TypeWrapper<char>( (char*)_store->address, _store->dims, _store->ndims, branchname  );

   if( Type == "l" ) 
      _ulong64_data[ branchname ]  = 
            new TypeWrapper<ULong64_t>( (ULong64_t*)_store->address, _store->dims, _store->ndims, branchname  );

   if( Type == "L" ) 
      _long64_data[ branchname ]  = 
            new TypeWrapper<Long64_t>( (Long64_t*)_store->address, _store->dims, _store->ndims, branchname  );

   // load the memory_map structure for future use
   _branch_space[ branchname ] = _store;
   _ttree->SetBranchAddress( branchname.c_str(), _store->address );


   delete _str;
}


void Builder::SetAllBranches( int verbose )
{
   SetVerbosity ( verbose );
   if( kBranchesSet ) return ;

   TTree * ltree;
   TFriendElement * fe;

/*
   _ttree->LoadTree(0);
   if( ListOfFriends )
   {
         fe = (TFriendElement*) ListOfFriends->At(0);
         ltree = fe->GetTree();
         ltree->LoadTree(0);
   }
*/

   SetAllBranches( _ttree, verbose );

    
   if( ListOfFriends )
   {
      std::cout << "Builder::SetAllBranches Adding branches from "
             << ListOfFriends->GetSize()
             << " friends." << std::endl;
      for( int i = 0 ; i < ListOfFriends->GetSize() ; i++ )
      {
         fe = (TFriendElement*) ListOfFriends->At(i);
         ltree = fe->GetTree();
         SetAllBranches( ltree, verbose );
      }
   }

   kBranchesSet = true ;

}

void Builder::SetAllBranches( TTree * ltree, int verbose )
{
   TObjArray * branches = ltree->GetListOfBranches();

   if( branches == NULL )
   {
      std::cout << "Error in Builder::SetAllBranches, No branches for tree " << ltree->GetName() << std::endl;
      return;
   }
   std::string name;

   int nBranches =  branches->GetEntries();
   
   for( int i = 0 ; i < nBranches ; i++ )
   {
       name =  (*branches)[i]->GetName();

       if( kVerbosity >= 1 ) std::cout << "Builder::SetAllBranches Setting branch >" << name  <<"<" << std::endl; 
       SetBranch(  name.c_str() );

       // restore branches - ugly but needed after calls to ltree->GetMaximum() 
       // in SetBranch
       branches = ltree->GetListOfBranches();

       if( kVerbosity >= 1 ) 
            std::cout << "Builder::SetAllBranches  .... set?    " << _branch_space.count( name )   << std::endl
                      <<   " size:  " << _branch_space[ name ]->size    << std::endl
                      <<   " ndims: " << _branch_space[ name ]->ndims   << std::endl
                      <<   " dims0: " << _branch_space[ name ]->dims[0] << std::endl
                      <<   " title: " << _branch_space[ name ]->title   << std::endl
                      <<   " Type:  " << _branch_space[ name ]->Type    << std::endl
                      <<" "<<  std::endl; 
   }  
   
}

std::vector< TString > Builder::GetListOfKeys()
{
   std::vector< TString > Keys;
   std::map< std::string, struct memory_map * >::iterator i; 

   for ( i = _branch_space.begin() ; i != _branch_space.end() ; i++ )
       Keys.push_back( i->first );

   return Keys;
} 

std::vector< TString > Builder::GetListOfBranchTitles()
{
   std::vector< TString > Titles;
   std::map< std::string, struct memory_map * >::iterator i; 

   for ( i = _branch_space.begin() ; i != _branch_space.end() ; i++ )
       Titles.push_back( i->second->title.c_str() );

   return Titles;
}



int Builder::GetSize( const char * branch )
{

   std::string branchname( branch );
 
   if( _branch_space.count( branchname ) ) 
      return _branch_space[ branchname ]->size ;
   
   return -1;
}


void Builder::Get( const TString key, TypeWrapper<int>& wrapper )
{
    if( _int_data.count(key) )
    {
       wrapper = *_int_data[key]; 
       if( !use_list.count( key ) )
           use_list[key] = 1 ; 
    }
    else 
      std::cout << " Warning Builder::Get " << key << " has not been registered. "
                << " There is no access to this variable " << std::endl; 
}


void Builder::Get( const TString key, TypeWrapper<unsigned>& wrapper )
{
    if( _unsigned_data.count(key) )
    {
       wrapper = *_unsigned_data[key]; 

       if( !use_list.count( key ) )
           use_list[key] = 1 ; 
    }
    else 
      std::cout << " Warning Builder::Get " << key << " has not been registered. "
                << " There is no access to this variable " << std::endl; 
}

void Builder::Get( const TString key, TypeWrapper<float>& wrapper )
{

    if( _float_data.count( key ) )
    {
       wrapper =  *_float_data[key]; 
       if( !use_list.count( key ) )
           use_list[key] = 1 ; 
    }
    else 
      std::cout << " Warning Builder::Get " << key << " has not been registered. "
                << " There is no access to this variable " << std::endl; 
}

void Builder::Get( const TString key, TypeWrapper<double>& wrapper )
{
    if( _double_data.count( key ) )
    {

       wrapper =  *_double_data[key]; 

       if( !use_list.count( key ) )
           use_list[key] = 1 ; 
    }
    else 
      std::cout << " Warning Builder::Get " << key << " has not been registered. "
                << " There is no access to this variable " << std::endl; 
}



//
//  Be Careful, this routine 
//  breaks the friendship relationship within ROOT
//  if you are friends with a tree(or chain) that has more entries 
//  than the first tree in the chain.
//
//  This situation arises if you make a TChain from a series of 
//  of TTrees but make a single tree (designed for friendship later)
//  from the entire chain:
//  orig. vars: |---tree 1 ---|---tree2 -|---tree3 ----|
//  new vars:   |------------friend chain--------------|
//
//  I anticipate this happening when, say, you want to add 
//  upmu variables to the standard existing atmpd ntuples
//  but only want one output file.
//
//  Its not clear if this situation should be avoided, but 
//  I think its better than dragging around 500 core files
//  plus 500X files, where is the number of additional 
//  variable sets you have
//
void Builder::GetEntry( unsigned int x )
{

   // in the case above I think it is not possible
   // to avoid a double call to GetEntry in the friend 
   // based on the guts of ROOT.
   // I need to look closer as this is _stupid_.
   _ttree->GetEntry(x);
   fCurrentEntry = x ;

   TTree * ltree;
   TFriendElement * fe;
    
/*
   if( ListOfFriends && _ttree->IsA() == TChain::Class() )
   {
      for( unsigned int i = 0 ; i < ListOfFriends->GetSize() ; i++ )
      if( LongFriendList[i] == 1 ) 
      {
         fe = (TFriendElement*) ListOfFriends->At(i);
         ltree = fe->GetTree();
         ltree->GetEntry(x);
      }
   }

*/


}



void Builder::EnableUsedListOnly( )
{

  std::map< TString , int >::iterator i;

  // disable all branches
  _ttree->SetBranchStatus("*",0); 

  // renable the ones we have "gotten" above
  for( i = use_list.begin() ; i != use_list.end() ; i++ )
    _ttree->SetBranchStatus( i->first.Data(), 1 );

}


void Builder::AlignmentTest( AlignFriend * af , unsigned skip )
{

  if( !ListOfFriends )
  {
     std::cout << " Builder::AlignmentTest  No friends specified for alignment test." << std::endl;
     return ;
  }
  std::cout << " Builder::AlignmentTest  Now testing friend alignment. " << std::endl;

  std::stringstream ss;
 
  unsigned size = ListOfFriends->GetSize();

  ULong64_t parent;
  std::vector< ULong64_t * >  FriendIDs;
  std::vector< TString   * >  TreeNames;  


  TTree * ltree;
  TFriendElement * fe;
  _ttree->LoadTree(0);

  unsigned i = 0;
  for( i = 0 ; i < size ; i++ )
  {
     fe = (TFriendElement*) ListOfFriends->At(i);
     ltree = fe->GetTree();

     // this is the default for alignment IDs
     ss.str(""); ss << ltree->GetName() << "_align" ;

     if ( ! Existence( ss.str().c_str() ) )
     {
        std::cout << "  warning, no aligmnent ID found in " << ltree->GetName() << std::endl; 
        std::cout << "    alignment check on this friend tree will be skipped " << std::endl;
        continue ;
     }

     std::cout << "   Adding: " << ss.str().c_str() << " to list of variables for alignment check. " << std::endl;
     FriendIDs.push_back(  (ULong64_t*) this->Get( ss.str().c_str()  ) );
     TreeNames.push_back(  new TString ( ltree->GetName() ) ) ;
     
  }// end of check for friends


  unsigned nEntries = GetEntries() ;
  
  bool kAligned = true;
  for ( unsigned i = 0 ; i < nEntries ; i += skip  )
  {

    this->GetEntry(i);

    parent = af->BuildAlignmentVar();
//  std::cout << parent << " : "  ;
    for ( unsigned j = 0 ; j < FriendIDs.size()  ; j++ )
    {
       kAligned =  *FriendIDs[j] == parent;
//     std::cout << *FriendIDs[j] << " "  ;
       
       if( !kAligned )
       {	
         std::cout << "Misalignment at entry " << i << " from tree: " << TreeNames[i]->Data() << std::endl;  
         std::cout << "   parent ID: " << parent << "  friend ID:" << *FriendIDs[j] << std::endl; 
         std::cout << "   returning." << std::endl;
         return;
       } 
    
    }
//  std::cout << std::endl;

  }




}


TArrayF * Builder::GetTArrayF( const TString key )
{
    if( _tarrayf_data.count(key) )
    {
       if( !use_list.count( key ) )
           use_list[key] = 1 ; 

       return _tarrayf_data[key] ;
    }
    else 
      std::cout << " Warning Builder::GetTArrayF " << key << " has not been registered. "
                << " There is no access to this variable. Returning null pointer! " << std::endl; 
    return 0 ;
}



