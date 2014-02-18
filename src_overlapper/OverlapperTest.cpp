#include "iostream"
#include "stdio.h"
#include "string"
#include "vector"
#include "cstdlib"

#include <map>
#include <math.h>
#include "memory"
#include <algorithm>
#include <fstream>
#include "sstream"
#include "list"
#include "stdlib.h"
#include "time.h"
#include <stdint.h>
#include "IntervalFileOperations.h"
#include "ReadsOperation.h"
using namespace std;




int main(int argc, char* argv[])
{

	string IntervalFile;
	string reads_info_name="TransformedReads.txt";
	
	for(int i=1;i<argc;++i)
	{
		if(strcmp(argv[i],"Interval")==0)
		{
			i++;
			IntervalFile=(argv[i]);
			continue;
		}

		if(strcmp(argv[i],"TransformedReads")==0)
		{
			i++;
			reads_info_name=(argv[i]);
			continue;
		}

	}

	ifstream intervals_in(IntervalFile.c_str());
	string temp;

	vector<vector<int> > pb_sr_map;
	string tag,n_tag;
	

	/*
	ofstream out_transformed_reads("TransformedReads.txt");
	while(get_pb_intervals(intervals_in, tag, pb_sr_map, n_tag))
	{
	
		if(pb_sr_map.size()>1)
		{
			out_transformed_reads<<tag<<endl;
			out_transformed_reads<<pb_sr_map.size()<<endl;

			for(int j=0;j<pb_sr_map.size();++j)
			{
				out_transformed_reads<<pb_sr_map[j].size()<<endl;
				for(int k=0;k<pb_sr_map[j].size();++k)
				{
					out_transformed_reads<<pb_sr_map[j][k]<<" ";
				}
				out_transformed_reads<<endl;
			}
			
		}

		pb_sr_map.clear();
	}
	

	*/
	reads_table reads_table0;
	ConstructReadsOverlaps_pb( reads_info_name,&reads_table0);

}