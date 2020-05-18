#include <iostream>
#include <stdlib.h>
#include <map>
#include <vector>

//////////////////////
//
//   This code is designed to automagically create
// projected 2- and 1-D delta chi^2 distributions
// from root files created with the Osc3++ libraries
//
// This code will not compile without Osc3++/tools
// 
// axisTree and stepX2Tree are required
//  200911026 rvw
/////////////////////////


#include "TTree.h"
#include "TString.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TDirectory.h"

#include "DataManager.h"
#include <sstream>

int ErrorExit( const char * );


using namespace std;
int main(int argc, char * argv[] )
{
  if( argc < 2 )
  {
     std::cerr << std::endl;
     std::cerr << std::endl;
     std::cerr << "Usage: ./BuildContourPlots [infile.root] "<< std::endl;
     std::cerr << "       rootfile must contain axisTree and stepX2Tree! "<< std::endl;

     std::cerr << std::endl;
     std::cerr << std::endl;
     std::cerr << "      You may alternatively specify which MC pulls should be saved: " << std::endl; 
     std::cerr << "      ./BuildContourPlots [infile.root] [point] "<< std::endl;

     std::cerr << std::endl;
     std::cerr << std::endl;
     std::cerr << "      You may additionally specify MC generation info : " << std::endl; 
     std::cerr << "      ./BuildContourPlots [infile.root] [point] [genpoint] "<< std::endl;

     std::cerr << std::endl;
     std::cerr << std::endl;
     exit( -1 );
  }

  int nMax = 400;
  unsigned int i, j;
  unsigned int k; 

  double minSyst[nMax]; 
  int nBins;


  TFile *_file = new TFile(argv[1] );   // the result of several fit files added together

  int WhichPoint = ( argc < 3 ? -1 : atoi( argv[2] ) );
  int GenPoint   = ( argc < 4 ?  0 : atoi( argv[3] ) );

  std::string outfile = "ChiSquared.root" ; 
  if ( argc >= 5 )  outfile = argv[4];


  TTree * axisTree   = (TTree*) _file->Get("axisTree"   );  
  TTree * stepX2Tree = (TTree*) _file->Get("stepX2Tree" );  
  TTree * errorTree  = (TTree*) _file->Get("errorTree"  );  


  if( axisTree   == 0 ) return ErrorExit("axisTree"); 
  if( stepX2Tree == 0 ) return ErrorExit("stepX2Tree");
  if( errorTree  == 0 ) return ErrorExit("errorTree"); 


  DataManager * axes = new DataManager( axisTree   );
  DataManager * X2t  = new DataManager( stepX2Tree );

  axes->SetAllBranches( );
  X2t ->SetAllBranches( );

  TFile * Out = new TFile( outfile.c_str() ,"recreate");

  double   gX2;
  double   mX2         = 1.e25; // the minimum chi2
  double * lX2         = (double*) X2t->Get("X2");
  int      nX;
  int      nY;
  double * Xaxis;
  double * Yaxis;

  double * Systematic  = (double*) X2t ->Get("Epsilons");
  int    * Point       = (  int *) X2t ->Get("Point"); 
  int    * nErrors     = (  int *) axes->Get("nErrors");
  double * Sigmas      = (double*) axes->Get("Sigmas");
  double * mVar;
  int      mPoint;

  axisTree->GetEntry(0);  // load variables

  std::vector< TString > axesBranches = axes->GetListOfKeys();
  std::vector< TString > X2tBranches  = X2t ->GetListOfKeys();

  std::vector< TString > axesBranchTitles = axes->GetListOfBranchTitles();
  std::vector< TString > X2tBranchTitles  = X2t ->GetListOfBranchTitles();

  std::stringstream ss;   
  std::stringstream title;   


  std::vector< TH2D * > ListOf2D;   
  std::vector< TH1D * > ListOf1D;   
  std::map< TString, double > ListOfMins;   
  std::map< TString, double > ListOfGens;   

  std::map< TString, double * > ListOfArrayMins;   
  std::map< TString, double * > ListOfArrayGens;   

  std::map< unsigned long, std::vector<TH2D*> > ListOf2DChildren;
  std::map< unsigned long, std::vector<TH1D*> > ListOf1DChildren;
  std::map< unsigned long, std::vector<TH1D*> > ListOf1DErrorChildren;
 

  TString * ErrorLabel = new TString("dummy"); 
  errorTree->SetBranchAddress( "ErrorLabel" , &ErrorLabel );

  TH1D * __minSyst  = new TH1D("SystMin" ,"", *nErrors, 0, *nErrors); 
  for( j = 0 ; j <*nErrors; j++ )
  {
     errorTree->GetEntry(j);

     __minSyst->GetXaxis()->SetBinLabel( j+1 , ErrorLabel->Data() ); 
     __minSyst->SetBinContent( j+1 , minSyst[j] );

  } 
  __minSyst->SetTitle("Fitted Errors #epsilon_{i} / #sigma_{i}");
  __minSyst->GetYaxis()->SetTitle("#epsilon_{i} / #sigma_{i}");
  __minSyst->GetXaxis()->SetTitle("");

  


  TH2D * l2D;
  TH2D * l2DChild;
  TH1D * l1DChild;
  TH1D * l1D;

  TString * XaxisName = 0;
  TString * YaxisName = 0;
  TString * ZaxisName = 0;
  TString * ParentName = 0;

  // Build all of the 1- and 2-D histograms that we will later fill
  for( i = 0 ; i < axesBranches.size() ; i++ )
    for( j = 0 ; j < axesBranches.size() ; j++ )
    {
    
       if( !axesBranches[i].Contains("axs") || !axesBranches[j].Contains("axs")  ) 
          continue;

       if( axesBranches[i].Contains("Length") || axesBranches[j].Contains("Length")  ) 
          continue;

       nX = axes->GetSize( axesBranches[i].Data() ) / sizeof(double) ; 
       nY = axes->GetSize( axesBranches[j].Data() ) / sizeof(double) ; 

       Xaxis = ( double*) axes->Get( axesBranches[i].Data() );		//here??? MJIANG
       Yaxis = ( double*) axes->Get( axesBranches[j].Data() );

//     if (Xaxis[0] < 0.) Xaxis[0] = 0;

       if( nX == 3 || nY == 3 ) continue;  // shall we skip axes that have only one bin

       if( XaxisName ) delete XaxisName;              
       if( YaxisName ) delete YaxisName;              

       XaxisName =  new TString( axesBranches[i].Data() );
       YaxisName =  new TString( axesBranches[j].Data() );
 
       XaxisName->ReplaceAll("axs", 3, "", 0 );
       YaxisName->ReplaceAll("axs", 3, "", 0 );

       ss.str(""); ss << "X2" << XaxisName->Data() << YaxisName->Data();
       title.str(""); title<< XaxisName->Data() << ":::" << YaxisName->Data();
       if( i != j )
       {
          l2D = new TH2D(ss.str().c_str(), title.str().c_str() , nX-2, Xaxis, nY-2, Yaxis ); 
          ListOf2D.push_back( l2D );

          if( ParentName ) delete ParentName;
          ParentName = new TString( ss.str().c_str() );

          // add us some chilluns!
          for( k = 0 ; k < axesBranches.size() ; k++ )
            if( k != i && k != j )
            {  
                if( !axesBranches[k].Contains("axs") || !axesBranches[k].Contains("axs")  ) 
                   continue;

                if( axesBranches[k].Contains("Length") || axesBranches[k].Contains("Length")  ) 
                   continue;

               if( ZaxisName ) delete ZaxisName;   
               ZaxisName = new TString( axesBranches[k].Data() );             
               ZaxisName->ReplaceAll("axs", 3, "", 0 );

               ss.str(""); ss << "X2" << XaxisName->Data() << YaxisName->Data() << "_" << ZaxisName->Data();
               title.str(""); title<< XaxisName->Data() << ":::" << YaxisName->Data() << " using " << ZaxisName->Data();
               
               l2DChild = new TH2D(ss.str().c_str(), title.str().c_str() , nX-2, Xaxis, nY-2, Yaxis ); 
               l2DChild->GetZaxis()->SetTitle( ZaxisName->Data());

               ListOf2DChildren[ l2D->Hash() ].push_back( l2DChild );

            }

       }

       // build minimum X2 histograms 
       if( i == j )
       {
          ss.str(""); ss << "X2min" << XaxisName->Data();
          l1D = new TH1D(ss.str().c_str(),XaxisName->Data(), nX-2, Xaxis);
          ListOf1D.push_back( l1D );
          ListOfMins[*XaxisName] = 1.0e10;

          // add 1-D children
          for( k = 0 ; k < axesBranches.size() ; k++ )
            if( k != i ) 
            {  
                if( !axesBranches[k].Contains("axs") || !axesBranches[k].Contains("axs")  ) 
                   continue;

                if( axesBranches[k].Contains("Length") || axesBranches[k].Contains("Length")  ) 
                   continue;

               if( ZaxisName ) delete ZaxisName;   
               ZaxisName = new TString( axesBranches[k].Data() );             
               ZaxisName->ReplaceAll("axs", 3, "", 0 );

               ss.str(""); ss << "X2min" << XaxisName->Data() << "_" << ZaxisName->Data();
               title.str(""); title<< XaxisName->Data() << "   using " << ZaxisName->Data();
               
               l1DChild = new TH1D(ss.str().c_str(), title.str().c_str() , nX-2, Xaxis);
               l1DChild->GetYaxis()->SetTitle( ZaxisName->Data());

               ListOf1DChildren[ l1D->Hash() ].push_back( l1DChild );


            }
            else 
            {

               for( unsigned q = 0 ; q <*nErrors;  q++ )
               {
                  errorTree->GetEntry(q);

                  ss.str(""); ss << "X2min" << XaxisName->Data() << "_" << ErrorLabel->Data();
                  title.str(""); title<< XaxisName->Data() << "   using " << ErrorLabel->Data() ;

                  if ( ErrorLabel->Contains("error_number") )
                  {
                    ss <<  "_" << q ;
                    title<< "_" << q ;
                  }
               
                  l1DChild = new TH1D(ss.str().c_str(), title.str().c_str() , nX-2, Xaxis);

                  l1DChild->GetYaxis()->SetTitle( "#epsilon / #sigma " );
                  ListOf1DErrorChildren[ l1D->Hash() ].push_back( l1DChild );
               }
 
 
            }

       }// end of building X2 min histos

    }// end of loop over Key names


  int minPoint = -1;
  cout << "Pulls will be at " << WhichPoint << endl;
  for( i = 0 ; i < stepX2Tree->GetEntries() ; i++ ){
      
     stepX2Tree->GetEntry(i);

 
     if ( *Point == WhichPoint )
       for( j = 0 ; j <*nErrors; j++ )
          minSyst[j] = Systematic[j]/Sigmas[j];

     
 
     if ( *lX2 < mX2 && *lX2 > -1.0e-150 )  // found a new minimum!
     {	
       mX2    = * lX2   ;
       mPoint = * Point ;
      
       for( j = 0 ; j < X2tBranches.size() ; j++ )
       {
          mVar = (double*) X2t->Get( X2tBranches[j].Data() ); 
          if( X2tBranchTitles[j].Contains("[") ) 
            ListOfArrayMins[ X2tBranches[j] ] = mVar;
          else
            ListOfMins[ X2tBranches[j] ] = *mVar;          
       }

       // Whichpoint == -1 means to fill the 
       // systematic errors at the best fit point 
       if ( WhichPoint == -1 ) 
         for( j = 0 ; j <*nErrors; j++ )
            minSyst[j] = Systematic[j]/Sigmas[j];

        //std::cout << "Point " << i << " " << Systematic[0] << " "  << mX2 << std::endl;

     }// end of new minimum

     //for Feldamn Cousins and Sensitivity studies
     if( *Point == GenPoint )
     {
       for( j = 0 ; j < X2tBranches.size() ; j++ )
       {
          mVar = (double*) X2t->Get( X2tBranches[j].Data() ); 
          if( X2tBranchTitles[j].Contains("[") ) 
            ListOfArrayGens[ X2tBranches[j] ] = mVar;
          else
            ListOfGens[ X2tBranches[j] ] = *mVar;          

       }
        gX2 = *lX2;
     }
    

  }// end loop on stepX2Tree entries

  std::cout << "Minimum Info - " << " X2 " << mX2 << std::endl;
  for( j = 0 ; j < X2tBranches.size() ; j++ )
    std::cout << "      " << X2tBranches[j].Data() << " : " << ListOfMins[X2tBranches[j]] <<std::endl;
  std::cout << " at point " << mPoint << std::endl;
  std::cout << std::endl;

  Out->cd();
  TTree * best = new TTree("bestfit", ""); 
  
  // Reload information from the best fit point  
  stepX2Tree->GetEntry( mPoint );

  // Best fit Point
  for( j = 0 ; j < X2tBranches.size() ; j++ )
  {
     ss.str(""); ss << X2tBranchTitles[j].Data() << "/D" ;
     if( X2tBranchTitles[j].Contains( "[" ) ) 
       best->Branch( X2tBranches[j].Data(),  &ListOfArrayMins[X2tBranches[j]][0], ss.str().c_str() ); 
     else 
       best->Branch( X2tBranches[j].Data(),  &ListOfMins[X2tBranches[j]], ss.str().c_str() ); 

  }


  // Generated Point
  for( j = 0 ; j < X2tBranches.size() ; j++ )
  {
     std::stringstream gen;

     gen.str(""); gen << "g" << X2tBranches[j].Data();
      ss.str(""); ss  << "g" << X2tBranchTitles[j].Data() ;
     if( X2tBranchTitles[j].Contains( "[" ) ) 
       best->Branch( gen.str().c_str(),  &ListOfArrayGens[X2tBranches[j]][0], ss.str().c_str() ); 
     else 
       best->Branch( gen.str().c_str(),  &ListOfGens[X2tBranches[j]], ss.str().c_str() ); 
  }

  best->Branch("mX2", &mX2 , "mX2/D" );
  best->Branch("mPoint", &mPoint , "mPoint/I" );
  best->Branch("gX2", &gX2 , "gX2/D" );
  best->Branch("gPoint", &GenPoint, "gPoint/I" );


  best->Fill();



  double Large = 1.0e10;
  int   ProposedBin;


  int nBinsX;
  int nBinsY;
  for( k = 0 ; k < ListOf2D.size() ; k++ )
  {
    l2D = ListOf2D[k];
    nBinsX = l2D->GetNbinsX(); 
    nBinsY = l2D->GetNbinsY(); 
    
    for( i = 1; i<= l2D->GetNbinsX(); i++)
       for( j = 1; j<= l2D->GetNbinsY(); j++)
         l2D->SetCellContent(i,j, Large);
  } 



  for( k = 0 ; k < ListOf1D.size() ; k++ )
  {
    l1D = ListOf1D[k];
    nBinsX = l1D->GetNbinsX(); 
    
    for( i = 1; i<= l1D->GetNbinsX(); i++)
       l1D->SetBinContent(i,Large);
  } 


  double * lVarX;
  double * lVarY;
  double * lVarZ;
  int len;
  int loc;  

  for( i = 0 ; i < stepX2Tree->GetEntries() ; i++ )
  {
     stepX2Tree->GetEntry(i);
       
     // local Chi2 is nonzero means it was fit
     // if( *lX2 < 1.0e-20 ) continue;
   
     for( k = 0 ; k < ListOf2D.size() ; k++ )
     {
       if( XaxisName ) delete XaxisName;
       if( YaxisName ) delete YaxisName;
        
       l2D = ListOf2D[k];
       XaxisName = new TString( l2D->GetTitle() );
       YaxisName = new TString( l2D->GetTitle() );

       loc = XaxisName->Index(":::");
       len = XaxisName->Length();
       XaxisName->Replace(loc,len-loc,"",0); 
       YaxisName->Replace(0,loc+3,"", 0);
 
//     if( i == 0 )  
//       std::cout << "Local TH2D's axes: " << XaxisName->Data() << " " << YaxisName->Data() << std::endl;

       lVarX = (double*) X2t->Get( XaxisName->Data() );
       lVarY = (double*) X2t->Get( YaxisName->Data() );

       ProposedBin = l2D->FindBin( *lVarX, *lVarY );
       if( (*lX2 - mX2) < l2D->GetBinContent( ProposedBin ) )
       {
          l2D->SetBinContent( ProposedBin, *lX2 - mX2 );

          for( j = 0 ; j < ListOf2DChildren[ l2D->Hash() ].size() ; j++ )
          {
              l2DChild = ListOf2DChildren[ l2D->Hash() ][j];
              if( ZaxisName ) delete ZaxisName;
              ZaxisName = new TString( l2DChild->GetZaxis()->GetTitle() ); 

              lVarZ = (double*) X2t->Get( ZaxisName->Data() );
              l2DChild->SetBinContent( ProposedBin, *lVarZ );
          }


       }

     }// end of loop on 2d histograms
	    

     for( k = 0 ; k < ListOf1D.size() ; k++ )
     {
       if( XaxisName ) delete XaxisName;
        
       l1D = ListOf1D[k];
       XaxisName = new TString( l1D->GetTitle() );
 
       lVarX = (double*) X2t->Get( XaxisName->Data() );

       ProposedBin = l1D->FindBin( *lVarX  );
       if( (*lX2 - mX2) < l1D->GetBinContent( ProposedBin ) )
       {
          l1D->SetBinContent( ProposedBin, *lX2 - mX2 );

          for( j = 0 ; j < ListOf1DChildren[ l1D->Hash() ].size() ; j++ )
          {
              l1DChild = ListOf1DChildren[ l1D->Hash() ][j];
              if( ZaxisName ) delete ZaxisName;
              ZaxisName = new TString( l1DChild->GetYaxis()->GetTitle() ); 

              lVarZ = (double*) X2t->Get( ZaxisName->Data() );
              l1DChild->SetBinContent( ProposedBin, *lVarZ );
          }

          // Fill error histograms
          for( j = 0 ; j < ListOf1DErrorChildren[ l1D->Hash() ].size() ; j++ )
          {
              l1DChild = ListOf1DErrorChildren[ l1D->Hash() ][j];

              // find correct error. Error-holding histograms were added to the 
              // vector in the order [0,nerrors] so use this: 
              l1DChild->SetBinContent( ProposedBin, Systematic[j]/Sigmas[j] );
          }

       }
     }// end of loop on 1d histograms


  }// end of loop on entries of stepX2Tree 
             

 /// Write all of the output at once
 best->Write();

 // now clone and write the errorTree with no midifications to the output file
 TTree * errorTreeClone = (TTree*) errorTree->CloneTree();
 errorTreeClone->Write();

 for( k = 0 ; k < ListOf2D.size() ; k++ )
   ListOf2D[k]->Write();

 for( k = 0 ; k < ListOf1D.size() ; k++ )
   ListOf1D[k]->Write();

 std::map< unsigned long, std::vector<TH2D *> >::iterator __i;
 for( __i = ListOf2DChildren.begin() ; __i != ListOf2DChildren.end() ; __i++ )
    for( i = 0 ; i < (__i->second).size() ; i++ )
       (__i->second)[i]->Write();

 std::map< unsigned long, std::vector<TH1D *> >::iterator __j;
 for( __j = ListOf1DChildren.begin() ; __j != ListOf1DChildren.end() ; __j++ )
    for( i = 0 ; i < (__j->second).size() ; i++ )
       (__j->second)[i]->Write();

 __minSyst->Write();

 TDirectory * sub =  Out->mkdir("systematic_minima");
 sub->cd();

 for( __j = ListOf1DErrorChildren.begin() ; __j != ListOf1DErrorChildren.end() ; __j++ )
    for( i = 0 ; i < (__j->second).size() ; i++ )
       (__j->second)[i]->Write();


  Out->Close();

  return 0;

}


int ErrorExit( const char  * treename )
{
    std::cout << "Input file does not contain TTree \"" << treename << "\" and cannot be used with this program." << std::endl;
    std::cout << "  Please use output from an Osc3++ *Fit or *Sensitivity Process. " << std::endl;
    return -1;
}

