#include "TFile.h"
#include "TTree.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <sys/stat.h>
int main(int argc, char* argv[]){
  if (argc != 2 && argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <string>" << std::endl;
        return 1;
    }
  int number_of_ev_perfile = 100;

  std::string input = argv[1];
  std::string name_file = input.substr(0,input.find('.'));
  //std::cout<<name_file<<std::endl;
  mkdir(name_file.c_str(),0777);
  mkdir((name_file+"/input").c_str(),0777);
  mkdir((name_file+"/output").c_str(),0777);
  mkdir((name_file+"/log").c_str(),0777);
  if(argc==3){
  std::stringstream os;
  os<<argv[2];
  os>>number_of_ev_perfile;
  }
  TFile f(("Root_File/"+input).c_str(),"read");
  TTree* tree=dynamic_cast<TTree*>(f.Get("apv_raw"));
  unsigned int nevt = tree->GetEntries();
  std::cout <<" # entries "<<nevt<<std::endl;

  int nfiles=nevt/number_of_ev_perfile+1;
  std::cout <<" Number of files "<<nfiles<<std::endl;
  for (unsigned int ii=0;ii<nfiles;ii++){
    std::stringstream os;
    os<<"runID_"<<std::setw(4)<<std::setfill('0')<<ii<<".root";
    std::cout <<"File "<<ii+1<<" "<<os.str()<<std::endl;
    unsigned int fev=ii*number_of_ev_perfile;
    unsigned int lev=std::min(nevt,(ii+1)*number_of_ev_perfile)-1;
    std::cout <<"Event range : "<<fev<<"--"<<lev<<std::endl;
    std::string name_out(name_file+"/input/"+os.str());
    TFile outp(name_out.c_str(),"RECREATE");
    auto pippo=dynamic_cast<TTree*>(tree->CloneTree(0));
    for (unsigned int jj=fev;jj<=lev;jj++){
      tree->GetEntry(jj);
      pippo->Fill();
    }
   
    pippo->AutoSave();
    outp.Close();
  }
}

