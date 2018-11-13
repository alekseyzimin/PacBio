#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <vector>
#include "src_mega_reads/get_super_read_sizes_cmdline.hpp"
#include <cstdlib>
#include<map>

int32_t main (int argc, char *argv[]){
cmdline_parse args;
args.parse (argc, argv);

std::map<long,long> kunitigsizes;
std::ifstream myfile;
std::string line;
myfile.open(args.kunitig_lenghts_file_arg);
  if (myfile.is_open()){
      while ( std::getline (myfile,line) ){
      char*tmpstr=new char[line.length()+1];
      strcpy(tmpstr,line.c_str());
      char*rest;
      long kunumber=strtol(strtok_r(tmpstr," ",&rest),NULL,10);
      long kusize=strtol(strtok_r(NULL," ",&rest),NULL,10);
      kunitigsizes[kunumber]=kusize;
      }
  }
  myfile.close();


long minkunitigsize=10000000;
for (std::map<long,long>::iterator it=kunitigsizes.begin(); it!=kunitigsizes.end(); ++it)
    if(it ->second < minkunitigsize)
        minkunitigsize=it ->second;

minkunitigsize--;

myfile.open(args.super_reads_file_arg);
  if (myfile.is_open()){
      while ( std::getline (myfile,line) ){
      long sr_size=0;
      char*tmpstr=new char[line.length()+1];
      strcpy(tmpstr,line.c_str());
      char*rest;
      long kunumber=strtol(strtok_r(tmpstr,"_FR",&rest),NULL,10);
      sr_size=kunitigsizes[kunumber];
      while(1){
        char*token=strtok_r(NULL,"_FR",&rest);
        if(token==NULL) break;
        long kunumber=strtol(token,NULL,10);
        sr_size+=(kunitigsizes[kunumber]-minkunitigsize);
      }
      std::cout <<line.c_str()<<" "<<sr_size<<std::endl;
      }
  }
  myfile.close();
}

