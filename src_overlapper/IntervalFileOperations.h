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
#include "BasicDataStructure.h"
using namespace std;

bool get_pb_intervals_orig(ifstream & interval_in, string &tag, vector<vector<int> >& pb_sr_map, string & n_tag, contigs_info *contigs_info)
{
	
	string temp,str;
	if(!getline(interval_in,temp))
	{return 0;}
	if(temp[temp.size()-1]=='\n'||temp[temp.size()-1]=='\r')
	{temp.resize(temp.size()-1);}
	if(temp.size()==0)
	{return 0;}
	int overlap = -1000;
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
					overlap = atoi(int_str.c_str());
					int_str.clear();//dist
					continue;
				}

				if(str[i]=='R')
				{
					ctg=-ctg;
				}

				ctgs_in_a_read.push_back(ctg);
				int_str.clear();

				if (ctgs_in_a_read.size() >= 2)
				{
					int ctg1 = ctgs_in_a_read[ctgs_in_a_read.size() - 2];
					int ctg2 = ctgs_in_a_read[ctgs_in_a_read.size() - 1];
					if (ctg1 < 0)
					{
						contigs_info->contig_adjacency_left[-ctg1][-ctg2].cov = 1;
						contigs_info->contig_adjacency_left[-ctg1][-ctg2].dist_sum = overlap;
					}
					else
					{
						contigs_info->contig_adjacency_right[ctg1][ctg2].cov = 1;
						contigs_info->contig_adjacency_right[ctg1][ctg2].dist_sum = overlap;
					}

					if (ctg2 > 0)
					{
						contigs_info->contig_adjacency_left[ctg2][ctg1].cov = 1;
						contigs_info->contig_adjacency_left[ctg2][ctg1].dist_sum = overlap;
					}
					else
					{
						contigs_info->contig_adjacency_right[-ctg2][-ctg1].cov = 1;
						contigs_info->contig_adjacency_right[-ctg2][-ctg1].dist_sum = overlap;
					}


				}

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
						overlap = atoi(int_str.c_str());
						int_str.clear();//dist
						int_str.clear();
						continue;
					}
					if(str[i]=='R')
					{ctg=-ctg;}

					ctgs_in_a_read.push_back(ctg);
					int_str.clear();


					if (ctgs_in_a_read.size() >= 2)
					{
						int ctg1 = ctgs_in_a_read[ctgs_in_a_read.size() - 2];
						int ctg2 = ctgs_in_a_read[ctgs_in_a_read.size() - 1];
						if (ctg1 < 0)
						{
							contigs_info->contig_adjacency_left[-ctg1][-ctg2].cov = 1;
							contigs_info->contig_adjacency_left[-ctg1][-ctg2].dist_sum = overlap;
						}
						else
						{
							contigs_info->contig_adjacency_right[ctg1][ctg2].cov = 1;
							contigs_info->contig_adjacency_right[ctg1][ctg2].dist_sum = overlap;
						}

						if (ctg2 > 0)
						{
							contigs_info->contig_adjacency_left[ctg2][ctg1].cov = 1;
							contigs_info->contig_adjacency_left[ctg2][ctg1].dist_sum = overlap;
						}
						else
						{
							contigs_info->contig_adjacency_right[-ctg2][-ctg1].cov = 1;
							contigs_info->contig_adjacency_right[-ctg2][-ctg1].dist_sum = overlap;
						}


					}
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


bool get_pb_intervals(ifstream & interval_in, string &tag, vector<vector<int> >& pb_sr_map, string & n_tag, contigs_info *contigs_info)
{

	string temp, str;
	if (!getline(interval_in, temp))
	{
		return 0;
	}
	if (temp[temp.size() - 1] == '\n' || temp[temp.size() - 1] == '\r')
	{
		temp.resize(temp.size() - 1);
	}
	if (temp.size() == 0)
	{
		return 0;
	}
	int overlap = -1000;
	bool flip = 0;
	str.clear();
	vector<int> ctgs_in_a_read;
	//if (temp[temp.size() - 1] == ':' || temp[temp.size() - 2] == ':')
	if (temp[0] == '>')
	{
		tag = temp;
	}
	else
	{
		tag = n_tag;
		str = temp;
		int ctg;
		string int_str;

		for (int i = 0; i<str.size(); ++i)
		{
			if (str[i] == ' ')
			{
				
				break;
			}
			if (str[i] <= '9'&&str[i] >= '0')
			{
				int_str.push_back(str[i]);
			}
			else
			{
				if (str[i] == 'F' || str[i] == 'R')
				{
					ctg = atoi(int_str.c_str());
				}
				else
				{
					overlap = atoi(int_str.c_str());
					int_str.clear();//dist
					continue;
				}

				if (str[i] == 'R')
				{
					ctg = -ctg;
				}

				ctgs_in_a_read.push_back(ctg);
				int_str.clear();

				if (ctgs_in_a_read.size() >= 2)
				{
					int ctg1 = ctgs_in_a_read[ctgs_in_a_read.size() - 2];
					int ctg2 = ctgs_in_a_read[ctgs_in_a_read.size() - 1];
					if (ctg1 < 0)
					{
						contigs_info->contig_adjacency_left[-ctg1][-ctg2].cov = 1;
						contigs_info->contig_adjacency_left[-ctg1][-ctg2].dist_sum = overlap;
					}
					else
					{
						contigs_info->contig_adjacency_right[ctg1][ctg2].cov = 1;
						contigs_info->contig_adjacency_right[ctg1][ctg2].dist_sum = overlap;
					}

					if (ctg2 > 0)
					{
						contigs_info->contig_adjacency_left[ctg2][ctg1].cov = 1;
						contigs_info->contig_adjacency_left[ctg2][ctg1].dist_sum = overlap;
					}
					else
					{
						contigs_info->contig_adjacency_right[-ctg2][-ctg1].cov = 1;
						contigs_info->contig_adjacency_right[-ctg2][-ctg1].dist_sum = overlap;
					}


				}

			}
		}

		pb_sr_map.push_back(ctgs_in_a_read);
		ctgs_in_a_read.clear();
	}


	while (getline(interval_in, temp))
	{
		

		if (temp[temp.size() - 1] == '\n' || temp[temp.size() - 1] == '\r')
		{
			temp.resize(temp.size() - 1);
		}

//		if (temp.size()>0 && (temp[temp.size() - 1] == ':' || temp[temp.size() - 2] == ':'))
		if (temp[0]=='>')
		{
			n_tag = temp;
			return 1;
		}
		else
		{
			str = temp;
			int ctg;
			string int_str;
			for (int i = 0; i<str.size(); ++i)
			{
				if (str[i] == ' ')
				{
					
					break;
				}
				if (str[i] <= '9'&&str[i] >= '0')
				{
					int_str.push_back(str[i]);
				}
				else
				{
					if (str[i] == 'F' || str[i] == 'R')
					{
						ctg = atoi(int_str.c_str());
					}
					else
					{
						overlap = atoi(int_str.c_str());
						int_str.clear();//dist
						int_str.clear();
						continue;
					}
					if (str[i] == 'R')
					{
						ctg = -ctg;
					}

					ctgs_in_a_read.push_back(ctg);
					int_str.clear();


					if (ctgs_in_a_read.size() >= 2)
					{
						int ctg1 = ctgs_in_a_read[ctgs_in_a_read.size() - 2];
						int ctg2 = ctgs_in_a_read[ctgs_in_a_read.size() - 1];
						if (ctg1 < 0)
						{
							contigs_info->contig_adjacency_left[-ctg1][-ctg2].cov = 1;
							contigs_info->contig_adjacency_left[-ctg1][-ctg2].dist_sum = overlap;
						}
						else
						{
							contigs_info->contig_adjacency_right[ctg1][ctg2].cov = 1;
							contigs_info->contig_adjacency_right[ctg1][ctg2].dist_sum = overlap;
						}

						if (ctg2 > 0)
						{
							contigs_info->contig_adjacency_left[ctg2][ctg1].cov = 1;
							contigs_info->contig_adjacency_left[ctg2][ctg1].dist_sum = overlap;
						}
						else
						{
							contigs_info->contig_adjacency_right[-ctg2][-ctg1].cov = 1;
							contigs_info->contig_adjacency_right[-ctg2][-ctg1].dist_sum = overlap;
						}


					}
				}
			}

			
			pb_sr_map.push_back(ctgs_in_a_read);
			ctgs_in_a_read.clear();

		}

	}
	return 1;
}










#endif
