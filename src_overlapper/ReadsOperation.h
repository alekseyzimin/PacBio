#ifndef __READS_OPERATION_H
#define __READS_OPERATION_H

#include <iostream>
#include <string>
#include <string.h>
#include <stdint.h>
#include <vector>
#include <map>
#include <list>
#include <algorithm>
#include <fstream>
#include "time.h"
#include "BasicDataStructure.h"


using namespace std;


struct SuperRead_t
{
	int PathsCnt;
	uint64_t idx;
	int PathsLength;
	int WellID;
	char seq[10000];
	bool PathFound;
	uint16_t depth;
	int ListSize;
	string extension;

};

struct search_info
{
	int pos1,pos2;
	bool RightSearch;
	bool Flip_End;
	int offset;
};


struct BFS_path_info_v2
{
	int cov;
	int depth;
	int len;
	bool BothEndsUsed;
	uint64_t last_kmer;
	struct edge_node* last_bkt_edge;
};


struct reads_overlap_info
{
	vector< map<int32_t, vector<int32_t> > > left_overlaps,right_overlaps;
	vector<int> cov_vt;
	vector<bool> contained_vt,used_vt,used_vt_left,used_vt_right;

};




struct BFS_reads_info
{
	int cov;
	int depth;
	int len;
	int last_read;
	vector<int> edge;
};



bool isSimplePath_read(reads_overlap_info *reads_overlap_info,int current_read,map<int,struct BFS_reads_info > &Visited_Path , map<int, int > &stacked_nodes)
{
	//return 1;
	int last_read=current_read;
	int dep=Visited_Path[current_read].depth;
	int node=last_read;
	for(int l=dep;l>=2;--l)//l>2
	{
	
		node=abs(Visited_Path[node].last_read);
		if(node==0)
		{return 1;}
		
		if(abs(stacked_nodes[node])>2)
		{return 1;}// only backtrack to here so return 1.
		
			
		if(stacked_nodes[node]>0&&reads_overlap_info->left_overlaps[abs(node)].size()>1)
		{return 0;}
		if(stacked_nodes[node]<0&&reads_overlap_info->right_overlaps[abs(node)].size()>1)
		{return 0;}
	}
	return 1;
}

void BreakLinks_read( reads_overlap_info *reads_overlap_info, map<int, int > &stacked_nodes,int node1, int node2)
{
	node1=abs(node1);
	node2=abs(node2);
	if(stacked_nodes[node1]>0)
	{
		int temp=abs(node2);
		if(stacked_nodes[node2]<0)
		{
			temp=-temp;
		}
		if(reads_overlap_info->right_overlaps[abs(node1)].count(temp))
		{
			reads_overlap_info->right_overlaps[abs(node1)].erase(temp);
			if(stacked_nodes[node1]>1)
			{
				stacked_nodes[abs(node1)]--;
			}
		}

	}


	if(stacked_nodes[node1]<0)
	{
		int temp=abs(node2);
		if(stacked_nodes[node2]>0)
		{
			temp=-temp;
		}

		if(reads_overlap_info->left_overlaps[abs(node1)].count(temp))
		{
			reads_overlap_info->left_overlaps[abs(node1)].erase(temp);
			if(stacked_nodes[node1]<-1)
			{
				stacked_nodes[abs(node1)]++;		
			}
		}


	}

	
	if(stacked_nodes[node2]<0)
	{
		int temp=abs(node1);
		if(stacked_nodes[node1]>0)
		{
			temp=-temp;
		}
		if(reads_overlap_info->right_overlaps[abs(node2)].count(temp))
		{
			reads_overlap_info->right_overlaps[abs(node2)].erase(temp);
				
		}
	}

	if(stacked_nodes[node2]>0)
	{
		int temp=abs(node1);
		if(stacked_nodes[node1]<0)
		{
			temp=-temp;
		}
		if(reads_overlap_info->left_overlaps[abs(node2)].count(temp))
		{
			reads_overlap_info->left_overlaps[abs(node2)].erase(temp);
				
		}
	}

}

void BacktrackBubbleRemoval_read(reads_overlap_info *reads_overlap_info,int last_read,int beg_read,map<int,struct BFS_reads_info > & Visited_Path , map<int ,int > &stacked_nodes)
{
	beg_read=abs(beg_read);
	last_read=abs(last_read);	
	int current_read=last_read;
	int dep=Visited_Path[last_read].depth;
	for(int l=dep;l>1;--l)
	{
		int previous_read=current_read;

		current_read=Visited_Path[current_read].last_read;
		if(current_read==0)
		{return;}
		
		if(stacked_nodes[current_read]>=1||stacked_nodes[current_read]<=-1)
		{
			

			if(abs(stacked_nodes[current_read])>2)
			{
				BreakLinks_read(reads_overlap_info,stacked_nodes,current_read,previous_read);
				break;
			}
			//else
			
			BreakLinks_read(reads_overlap_info,stacked_nodes,current_read,previous_read);
			
			if(Visited_Path[current_read].last_read==NULL)
			{
				break;
			}

			int free_read=current_read;

			if(beg_read==free_read)
			{
				break;
			}

			if(free_read!=last_read)
			{

				stacked_nodes[free_read]=stacked_nodes[free_read]/abs(stacked_nodes[free_read]);
				//freebkt->kmer_info.removed=1;			
			}
		}
	}
}


void BFSearchBubbleRemoval_read(reads_overlap_info *reads_overlap_info,reads_table *reads_table,int beg_read,int max_depth,int max_dist)
{
	map<int,struct BFS_reads_info > Visited_Path;
	map<int, int > stacked_nodes;
	int max_stack=300;
	int DepthTh=max_depth;
	int LenTh=30;
	bool RIGHT=0;
	map<int, list<int> > dist_reads;//neighborset
	dist_reads[0].push_back(beg_read);
	int NBs=1;
	int dist_searched=0;

	int new_node=abs(beg_read);
	if(beg_read>0)
	{
		stacked_nodes[abs(beg_read)]=1;
	}
	else
	{
		stacked_nodes[abs(beg_read)]=-1;
	}
	
	Visited_Path[new_node].cov=0;
	Visited_Path[new_node].depth=1;
	Visited_Path[new_node].last_read=0;
	
	map<int , list<int> >::iterator NB_it=dist_reads.begin();

	while(1)
	{
		NB_it=dist_reads.begin();
		if(NB_it==dist_reads.end())
		{break;}
		if(NB_it->second.size()==0)
		{dist_reads.erase(NB_it->first);continue;}

		if(NBs>max_stack)
		{
			break;
		}
		new_node=NB_it->second.front();

		NB_it->second.pop_front();
		NBs--;
		if(NB_it->second.size()==0)
		{
			dist_reads.erase(NB_it->first);
		}
		if(new_node>0)
		{
			RIGHT=1;
		}
		else
		{
			RIGHT=0;
		}

		new_node=abs(new_node);
		//if(new_node==173)
		//cout<<new_node<<endl;
		if(Visited_Path[new_node].depth>DepthTh||Visited_Path[new_node].len>LenTh)
		{continue;}
		//if(new_node==239)
		//{cout<<"";}

		if(RIGHT)
		{
			int rb=reads_overlap_info->right_overlaps[new_node].size();			
			if(stacked_nodes[new_node]==1&&rb>0)
			{
				stacked_nodes[new_node]=1+rb;
			}
			if(rb==0&&reads_overlap_info->used_vt_right[abs(new_node)]==0)
			{
				stacked_nodes[new_node]=2;
				
				if(Visited_Path[abs(new_node)].depth>2)// let's define a tip to be non-extendable
				{
					continue;
				}
				if(abs(new_node)==abs(beg_read))
				{					
					continue;
				}
				//tip end reached so backtrack to the branching position.
				if(!isSimplePath_read(reads_overlap_info,new_node,Visited_Path, stacked_nodes))
				{
					continue;
				}

				BacktrackBubbleRemoval_read(reads_overlap_info,new_node,beg_read,Visited_Path,stacked_nodes);
				if(beg_read>0)
				{
					reads_overlap_info->used_vt_right[abs(beg_read)]=1;
				}
				else
				{
					reads_overlap_info->used_vt_left[abs(beg_read)]=1;
				}
				stacked_nodes[new_node]=1;
				continue;

			}

			map<int32_t, vector<int32_t> >::iterator tmp_it,tmp_it_n;
			for(tmp_it=reads_overlap_info->right_overlaps[new_node].begin();tmp_it!=reads_overlap_info->right_overlaps[new_node].end();)
			{
				tmp_it_n=tmp_it;
				tmp_it_n++;
				int next_read=tmp_it->first;
				// not in stack
				if(stacked_nodes[abs(next_read)]==0)
				{
//					Visited_Path[*ptr].cov=(int)(Visited_Path[new_node].cov+edge_ptr->edge_cov);
					Visited_Path[abs(next_read)].cov=(int)(Visited_Path[new_node].cov+reads_overlap_info->cov_vt[abs(next_read)]);
					Visited_Path[abs(next_read)].depth=(Visited_Path[abs(new_node)].depth+1);
					//int cum_len=(int)(Visited_Path[abs(next_read)].len+1);
					int cum_len=Visited_Path[abs(next_read)].depth;
					//Visited_Path[abs(next_read)].len=cum_len;//(Visited_Path[new_node].len+edge_ptr->len+1);
					Visited_Path[abs(next_read)].last_read=new_node;
					
					if(next_read<0)
					{
						stacked_nodes[abs(next_read)]=-1;	
					}
					else
					{
						stacked_nodes[abs(next_read)]=1;
					}
					dist_reads[cum_len].push_back(next_read);
					NBs++;

				}
				else
				{
					if((stacked_nodes[abs(next_read)]>0&&next_read>0)||(stacked_nodes[abs(next_read)]<0&&next_read<0))
					{
						
						if((Visited_Path[abs(new_node)].cov+reads_overlap_info->cov_vt[abs(next_read)]<=Visited_Path[abs(next_read)].cov))//||(BackCheckLoop(*ptr,new_node,Visited_Path)==1)//loop
						{
							//backtrack if the same direction is found
							
							if(!isSimplePath_read(reads_overlap_info,new_node,Visited_Path, stacked_nodes))
							{
								tmp_it=tmp_it_n;
								continue;
							}
							
							//backtrack the current path, common operation in this search
							if(stacked_nodes[new_node]>2)
							{
								BreakLinks_read(reads_overlap_info,stacked_nodes, new_node, next_read);
								tmp_it=tmp_it_n;
								continue;
							}
							else
							{
								if(stacked_nodes[new_node]<-2)
								{
									BreakLinks_read(reads_overlap_info,stacked_nodes, new_node, next_read);
									
									tmp_it=tmp_it_n;
									continue;
								}
								else
								{
									int free_node=(new_node);

									if(free_node==beg_read)
									{
										BreakLinks_read(reads_overlap_info,stacked_nodes, beg_read, new_node);
									}
									
									//free the node and edge.
									BreakLinks_read(reads_overlap_info,stacked_nodes,new_node,next_read);
									BacktrackBubbleRemoval_read(reads_overlap_info,new_node,beg_read,Visited_Path,stacked_nodes);
									tmp_it=tmp_it_n;
									continue;
			//
								}
							}

						}
						else
						{
							//backtrack the original path, rare operation in this search, can lead to errors
							
							//BacktrackBubbleRemoval(ht,merge_ht,*ptr,beg_bkt,beg_bkt,Visited_Path , stacked_nodes,K_size);
							BacktrackBubbleRemoval_read(reads_overlap_info,next_read,beg_read,Visited_Path,stacked_nodes);
							
							//marginal 
							Visited_Path[abs(next_read)].cov=(int)(Visited_Path[abs(new_node)].cov+reads_overlap_info->cov_vt[abs(next_read)]);
							Visited_Path[abs(next_read)].depth=Visited_Path[abs(new_node)].depth+1;
							Visited_Path[abs(next_read)].len=Visited_Path[abs(new_node)].depth;
							//marginal 

							Visited_Path[abs(next_read)].last_read=abs(new_node);
							tmp_it=tmp_it_n;
							continue;
						}

					}
					else
					{

						//don't do anything,since both strands are visited.

					}

				}
				tmp_it=tmp_it_n;
			}
			
		}
		else
		{

			int lb=reads_overlap_info->left_overlaps[new_node].size();			
			if(stacked_nodes[new_node]==-1&&lb>0)
			{
				stacked_nodes[new_node]=-1-lb;
			}
			if(lb==0&&reads_overlap_info->used_vt_left[abs(new_node)]==0)
			{
				stacked_nodes[new_node]=-2;
				if(Visited_Path[abs(new_node)].depth>2)// let's define a tip to be non-extendable
				{
					continue;
				}
				
				if(abs(new_node)==abs(beg_read))
				{					
					continue;
				}
				//tip end reached so backtrack to the branching position.
				if(!isSimplePath_read(reads_overlap_info,new_node,Visited_Path, stacked_nodes))
				{
					continue;
				}

				BacktrackBubbleRemoval_read(reads_overlap_info,new_node,beg_read,Visited_Path,stacked_nodes);
				if(beg_read>0)
				{
					reads_overlap_info->used_vt_right[abs(beg_read)]=1;
				}
				else
				{
					reads_overlap_info->used_vt_left[abs(beg_read)]=1;
				}
				stacked_nodes[new_node]=-1;
				continue;

			}

			map<int32_t, vector<int32_t> >::iterator tmp_it,tmp_it_n;
			for(tmp_it=reads_overlap_info->left_overlaps[new_node].begin();tmp_it!=reads_overlap_info->left_overlaps[new_node].end();)
			{
				tmp_it_n=tmp_it;
				tmp_it_n++;
				int next_read=tmp_it->first;
				// not in stack
				if(stacked_nodes[abs(next_read)]==0)
				{
//					Visited_Path[*ptr].cov=(int)(Visited_Path[new_node].cov+edge_ptr->edge_cov);
					Visited_Path[abs(next_read)].cov=(int)(Visited_Path[new_node].cov+reads_overlap_info->cov_vt[abs(next_read)]);
					Visited_Path[abs(next_read)].depth=(Visited_Path[abs(new_node)].depth+1);
					//int cum_len=(int)(Visited_Path[abs(next_read)].len+1);
					int cum_len=Visited_Path[abs(next_read)].depth;
					//Visited_Path[abs(next_read)].len=cum_len;//(Visited_Path[new_node].len+edge_ptr->len+1);
					Visited_Path[abs(next_read)].last_read=new_node;
					
					if(next_read<0)
					{
						stacked_nodes[abs(next_read)]=1;	
					}
					else
					{
						stacked_nodes[abs(next_read)]=-1;
					}
					dist_reads[cum_len].push_back(-next_read);
					NBs++;

				}
				else
				{
					if((stacked_nodes[abs(next_read)]<0&&next_read>0)||(stacked_nodes[abs(next_read)]>0&&next_read<0))
					{
						
						if((Visited_Path[abs(new_node)].cov+reads_overlap_info->cov_vt[abs(next_read)]<=Visited_Path[abs(next_read)].cov))//||(BackCheckLoop(*ptr,new_node,Visited_Path)==1)//loop
						{

							//backtrack if the same direction is found

							if(!isSimplePath_read(reads_overlap_info,new_node,Visited_Path, stacked_nodes))
							{
								tmp_it=tmp_it_n;
								continue;
							}
							
							//backtrack the current path, common operation in this search
							if(stacked_nodes[new_node]>2)
							{
								BreakLinks_read(reads_overlap_info,stacked_nodes, new_node, next_read);
								tmp_it=tmp_it_n;
								continue;
							}
							else
							{
								if(stacked_nodes[new_node]<-2)
								{
									BreakLinks_read(reads_overlap_info,stacked_nodes, new_node, next_read);
									tmp_it=tmp_it_n;
									continue;
								}
								else
								{
									int free_node=(new_node);

									if(free_node==beg_read)
									{
										BreakLinks_read(reads_overlap_info,stacked_nodes, beg_read, new_node);
									}
									
									//free the node and edge.
									BreakLinks_read(reads_overlap_info,stacked_nodes,new_node,next_read);
									BacktrackBubbleRemoval_read(reads_overlap_info,new_node,beg_read,Visited_Path,stacked_nodes);
									tmp_it=tmp_it_n;
									continue;
								}
							}
						}
						else
						{
							//backtrack the original path, rare operation in this search, can lead to errors
							BacktrackBubbleRemoval_read(reads_overlap_info,next_read,beg_read,Visited_Path,stacked_nodes);
							
							//marginal 
							Visited_Path[abs(next_read)].cov=(int)(Visited_Path[abs(new_node)].cov+reads_overlap_info->cov_vt[abs(next_read)]);
							Visited_Path[abs(next_read)].depth=Visited_Path[abs(new_node)].depth+1;
							Visited_Path[abs(next_read)].len=Visited_Path[abs(new_node)].depth;
							//marginal 

							Visited_Path[abs(next_read)].last_read=abs(new_node);
							tmp_it=tmp_it_n;
							continue;
						}

					}
					else
					{

						//don't do anything,since both strands are visited.

					}


				}

				
				tmp_it=tmp_it_n;




			}

		}

	}

}



void ConstructReadsOverlaps_pb(string reads_info_name,reads_table* reads_table)
{
	// 
	reads_overlap_info reads_overlap_info;
	bool EXACT=1;
	ifstream in_reads_info;
	in_reads_info.open(reads_info_name.c_str());
	ofstream out_extended_reads("Extended_Sr.txt");
	string tag;
	int32_t num_reads=0,n_selected_reads=0;
	int n_ctgs;
	int n_sr_lines;
	vector<vector<int32_t> > read_contig_index;
	map<int32_t,vector<int32_t> > contig_in_reads;
	map<vector<int>,int> selected_reads;

	reads_overlap_info.cov_vt.push_back(0);

	while(getline(in_reads_info,tag))
	{
		
		n_selected_reads=0;
		selected_reads.clear();
		read_contig_index.clear();
		contig_in_reads.clear();
		reads_overlap_info.cov_vt.clear();
		reads_overlap_info.left_overlaps.clear();
		reads_overlap_info.right_overlaps.clear();
		reads_overlap_info.contained_vt.clear();
		reads_overlap_info.cov_vt.push_back(0);
		reads_overlap_info.used_vt.clear();
		reads_overlap_info.used_vt_left.clear();
		reads_overlap_info.used_vt_right.clear();

		if(tag.size()<=1)
		{
			continue;
		}
		out_extended_reads<<tag<<endl;
		num_reads++;
		in_reads_info>>n_sr_lines;
		for (int sr_line=0;sr_line<n_sr_lines;++sr_line)
		{
			in_reads_info>>n_ctgs;
			vector<int> contigs_vt,contigs_vt_rc;
			int ctg;
			// get the contigs into a vector
			for (int i=0;i<n_ctgs;++i)
			{
				in_reads_info>>ctg;
				contigs_vt.push_back(ctg);
			}
		
			contigs_vt_rc=contigs_vt;
			reverse(contigs_vt_rc.begin(),contigs_vt_rc.end());
		
			// compare and take the smaller vector
			for (int i=0;i<n_ctgs;++i)
			{
				contigs_vt_rc[i]=-contigs_vt_rc[i];
			}
			bool take_rc=0;	
			for (int i=0;i<n_ctgs;++i)
			{
				if(contigs_vt_rc[i]>contigs_vt[i])
				{
					take_rc=0;
					break;
				}
				if(contigs_vt_rc[i]<contigs_vt[i])
				{
					take_rc=1;
					break;
				}
			}
			if(take_rc)
			{
				contigs_vt=contigs_vt_rc;
			}
			if(selected_reads[contigs_vt]!=0)
			{
			
				reads_overlap_info.cov_vt[abs(selected_reads[contigs_vt])]++;
				continue;
			}
			else
			{
				n_selected_reads++;
				if(take_rc)
				{
					selected_reads[contigs_vt]=-n_selected_reads;
					reads_overlap_info.cov_vt.push_back(1);
				}
				else
				{
					selected_reads[contigs_vt]=n_selected_reads;
					reads_overlap_info.cov_vt.push_back(1);
				}
			}

			//cout<<n_selected_reads<<endl;
			// now in selected reads
			for (int i=0;i<n_ctgs;++i)
			{
				ctg=contigs_vt[i];	
				if(contig_in_reads[abs(ctg)].size()==0||contig_in_reads[abs(ctg)].back()!=n_selected_reads)//check and don't save duplicate values
				{
					contig_in_reads[abs(ctg)].push_back(n_selected_reads);
				}
			}

			read_contig_index.push_back(contigs_vt);
	
		}

		///////for each pb read



		reads_overlap_info.left_overlaps.resize(n_selected_reads+1);
		reads_overlap_info.right_overlaps.resize(n_selected_reads+1);
		reads_overlap_info.contained_vt.resize(n_selected_reads+1);
		reads_overlap_info.used_vt.resize(n_selected_reads+1);
		reads_overlap_info.used_vt_left.resize(n_selected_reads+1);
		reads_overlap_info.used_vt_right.resize(n_selected_reads+1);
		
		for (int i=0;i<read_contig_index.size();++i)
		{

			if(EXACT)
			{

			
				for (int round1=1;round1<=2;++round1)
				{
					//front/rear overlap
					vector<int32_t> current_read=read_contig_index[i];
				
					if(round1==2)
					{
						reverse(current_read.begin(),current_read.end());
						for (int j=0;j<current_read.size();++j)
						{
							current_read[j]=-current_read[j];
						}

					}

			
					for (int r=0;r<contig_in_reads[abs(current_read[0])].size();++r)
					{
						int read_idx=contig_in_reads[abs(current_read[0])][r];
						if(read_idx-1==i)
						{
							continue;
						}

						vector<int32_t> ctgs_vt=read_contig_index[read_idx-1];
						for (int round2=1;round2<=2;++round2)
						{
							if(round2==2)
							{
								reverse(ctgs_vt.begin(),ctgs_vt.end());
								for (int j=0;j<ctgs_vt.size();++j)
								{
									ctgs_vt[j]=-ctgs_vt[j];
								}
							}
							for (int j=0;j<ctgs_vt.size();++j)
							{
								int pos=0;
								bool overlap=1;
								for(int k=j;k!=ctgs_vt.size();++k)
								{
									if (current_read[pos]!=ctgs_vt[k])
									{
										overlap=0;
										break;
									}
									pos++;
									if(pos==current_read.size()&&k!=ctgs_vt.size()-1)
									{
										overlap=0;
										break;
									}
								}
								//check if we have a proper overlap.
								if(overlap==1)
								{
									vector<int> edge;
									if(ctgs_vt.size()-j==current_read.size())
									{
										reads_overlap_info.contained_vt[i+1]=1;
										cout<<"";//contained overlap
									}


									for(int k=0;k<j;++k)
									{
										edge.push_back(ctgs_vt[k]);
									}

									// record the overlap and break;
									if(round1==1)
									{
										//front overlap
										if(round2==1)
										{
											//
											reads_overlap_info.left_overlaps[i+1][read_idx]=edge;
										}
										else
										{
											//flip
											
											reads_overlap_info.left_overlaps[i+1][-read_idx]=edge;
										}

									}
									else
									{
										//rear overlap

										reverse(edge.begin(),edge.end());
										for(int k=0;k<edge.size();++k)
										{
											edge[k]=-edge[k];
										}

										if(round2==1)
										{
											//flip
											
											reads_overlap_info.right_overlaps[i+1][-read_idx]=edge;
										}
										else
										{
										
											reads_overlap_info.right_overlaps[i+1][read_idx]=edge;
										}
									}

									break;

								}
					
							}
					
						}

					}

				}

			}

		}
	
		for(int i=0;i<reads_overlap_info.right_overlaps.size();++i)
		{
			map<int32_t, vector<int32_t> >::iterator temp_it1,temp_it2;
			for(temp_it1=reads_overlap_info.right_overlaps[i].begin();temp_it1!=reads_overlap_info.right_overlaps[i].end(); )
			{
				temp_it2=temp_it1;
				temp_it2++;
				if(reads_overlap_info.contained_vt[abs(temp_it1->first)]||reads_overlap_info.contained_vt[i])
				{
					reads_overlap_info.right_overlaps[i].erase(temp_it1);
				}
				temp_it1=temp_it2;
			}

			for(temp_it1=reads_overlap_info.left_overlaps[i].begin();temp_it1!=reads_overlap_info.left_overlaps[i].end(); )
			{
				temp_it2=temp_it1;
				temp_it2++;
				if(reads_overlap_info.contained_vt[abs(temp_it1->first)]||reads_overlap_info.contained_vt[i])
				{
					reads_overlap_info.left_overlaps[i].erase(temp_it1);
				}
				temp_it1=temp_it2;
			}

		}

		

		int max_dist=20;
		int max_depth=20;
	
		for(int i=1;i<reads_overlap_info.contained_vt.size();++i)
		{
			//cout<<i<<endl;
			//cout<<reads_overlap_info.right_overlaps[11].size()<<endl;
			int beg_read=i;
			

			if(reads_overlap_info.contained_vt[i]==0)
			{
				if(reads_overlap_info.right_overlaps[i].size()>1)
				{
					BFSearchBubbleRemoval_read(&reads_overlap_info,reads_table, beg_read, max_depth, max_dist);
				}
				//cout<<"l"<<endl;
				if(reads_overlap_info.left_overlaps[i].size()>1)
				{
					BFSearchBubbleRemoval_read(&reads_overlap_info,reads_table, -beg_read, max_depth, max_dist);
				}
			}
		}


		

		for(int read_idx=1;read_idx<reads_overlap_info.contained_vt.size();++read_idx)
		{
			vector<int> left_ext,right_ext,extended_read;
			if(reads_overlap_info.contained_vt[read_idx]==0&&reads_overlap_info.used_vt[read_idx]==0)
			{
				

				for (int it=1;it<=2;++it)
				{
					int current_read=read_idx;
					
					reads_overlap_info.used_vt[current_read]=1;
					bool Right;
					if(it==1)
					{
						Right=1;
					}
					else
					{
						Right=0;
					}

					while(1)
					{
						int next_read;
						if(Right==1)
						{		
							if(reads_overlap_info.right_overlaps[abs(current_read)].size()==1)
							{
								next_read=reads_overlap_info.right_overlaps[abs(current_read)].begin()->first;
								vector<int> edge=reads_overlap_info.right_overlaps[abs(current_read)].begin()->second;
								for(int e=0;e<edge.size();++e)
								{
									if(it==1)
									{
										right_ext.push_back(edge[e]);
									}
									else
									{
										left_ext.push_back(-edge[e]);
									}
								}
							}
							else
							{break;}
						}
						else
						{
							if(Right==0)
							{
								if(reads_overlap_info.left_overlaps[abs(current_read)].size()==1)
								{
									next_read=reads_overlap_info.left_overlaps[abs(current_read)].begin()->first;
									vector<int> edge=reads_overlap_info.left_overlaps[abs(current_read)].begin()->second;
									
									for(int e=0;e<edge.size();++e)
									{
										if(it==1)
										{
										right_ext.push_back(-edge[e]);
										}
										else
										{
											left_ext.push_back(edge[e]);
										}
									}
								}
								else
								{
									break;
								}
							}
						}
						
						if(next_read<0)
						{
							Right=!Right;
						}
						current_read=abs(next_read);
						if(reads_overlap_info.used_vt[abs(current_read)])
						{break;}
						reads_overlap_info.used_vt[abs(current_read)]=1;
					}
				}

				if(left_ext.size()>0||right_ext.size()>0)
				{
					for(int jj=0;jj<left_ext.size();++jj)
					{
						extended_read.push_back(left_ext[jj]);
						out_extended_reads<<left_ext[jj]<<" ";
					}
					for(int jj=0;jj<read_contig_index[read_idx-1].size();++jj)
					{
						extended_read.push_back(read_contig_index[read_idx-1][jj]);
						out_extended_reads<<read_contig_index[read_idx-1][jj]<<" ";
					}
					for(int jj=0;jj<right_ext.size();++jj)
					{
						extended_read.push_back(right_ext[jj]);
						out_extended_reads<<right_ext[jj]<<" ";
					}
					out_extended_reads<<endl;
				}

			}

			

		}
		cout<<"";
		/////////////finished




	}
	in_reads_info.close();
	out_extended_reads.close();

	cout<<"Overlap graph constructed"<<endl;
}






#endif
