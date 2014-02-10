#ifndef __INTERVAL_FILE_OPERATIONS_H
#define __INTERVAL_FILE_OPERATIONS_H




#include <iostream>
#include <string>
#include <string.h>
#include <stdint.h>
#include <vector>
#include <map>
#include <list>
#include <algorithm>
#include <fstream>
#include <sstream>
#include "time.h"
using namespace std;

bool get_pb_intervals(ifstream & interval_in, string &tag, vector<vector<int> >& pb_sr_map, string & n_tag)
{
	
	string temp,str;
	if(!getline(interval_in,temp))
	{return 0;}
	if(temp[temp.size()-1]=='\n'||temp[temp.size()-1]=='\r')
	{temp.resize(temp.size()-1);}
	if(temp.size()==0)
	{return 0;}

	bool flip=0;
	str.clear();
	vector<int> ctgs_in_a_read;
	if(temp[temp.size()-1]==':'||temp[temp.size()-2]==':')
	{tag=temp;}
	else
	{
		tag=n_tag;
		str=temp;
		int ctg;
		string int_str;
		for (int i=1;i<str.size();++i)
		{
			if(str[i]==' ')
			{
				if(str[i+2]=='+')
				{
					flip=0;
				}
				else
				{
					flip=1;
				}
				break;
			}
			if(str[i]<='9'&&str[i]>='0')
			{
				int_str.push_back(str[i]);
			}
			else
			{
				if(str[i]=='F'||str[i]=='R')
				{ctg=atoi(int_str.c_str());}
				else
				{
					int_str.clear();
					continue;
				}
				if(str[i]=='R')
				{ctg=-ctg;}

				ctgs_in_a_read.push_back(ctg);
				int_str.clear();
			}
		}

		if(flip==1)
		{
			reverse(ctgs_in_a_read.begin(),ctgs_in_a_read.end());
			for(int j=0;j<ctgs_in_a_read.size();++j)
			{
				ctgs_in_a_read[j]=-ctgs_in_a_read[j];
			}
		}

		pb_sr_map.push_back(ctgs_in_a_read);
		ctgs_in_a_read.clear();
	}

	
	while(getline(interval_in,temp))
	{
		if(temp[0]!='>')
		{continue;}
		
		if(temp[temp.size()-1]=='\n'||temp[temp.size()-1]=='\r')
		{temp.resize(temp.size()-1);}

		if(temp.size()>0&&(temp[temp.size()-1]==':'||temp[temp.size()-2]==':'))
		{
			n_tag=temp;
			return 1;
		}
		else
		{
			str=temp;
			int ctg;
			string int_str;
			for (int i=1;i<str.size();++i)
			{
				if(str[i]==' ')
				{
					if(str[i+2]=='+')
					{

						flip=0;
					}
					else
					{
						flip=1;
					} 
				

					break;
				}
				if(str[i]<='9'&&str[i]>='0')
				{
					int_str.push_back(str[i]);
				}
				else
				{
					if(str[i]=='F'||str[i]=='R')
					{ctg=atoi(int_str.c_str());}
					else
					{
						int_str.clear();
						continue;
					}
					if(str[i]=='R')
					{ctg=-ctg;}

					ctgs_in_a_read.push_back(ctg);
					int_str.clear();
				}
			}

			if(flip==1)
			{
				reverse(ctgs_in_a_read.begin(),ctgs_in_a_read.end());
				for(int j=0;j<ctgs_in_a_read.size();++j)
				{
					ctgs_in_a_read[j]=-ctgs_in_a_read[j];
				}
			}
			pb_sr_map.push_back(ctgs_in_a_read);
			ctgs_in_a_read.clear();
			
		}
		
	}
	return 1;
}












#endif
