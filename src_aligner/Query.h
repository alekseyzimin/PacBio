#ifndef __QUERY_H
#define __QUERY_H
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
#include "Alignment.h"
#include "BasicDataStructure.h"
#include "GraphConstruction.h"

using namespace std;


void Sparse_Kmer_Ref_Graph_Query(struct ref_t *ref, struct read_t *read,struct hashtable *ht,int K_size,int gap, struct qry_ret *qry_ret)//,map<uint64_t,int > &key_map,int &map_match)
{
	qry_ret->matched_seq.clear();
	qry_ret->match_cnt.clear();
	qry_ret->matched_seq.clear();
	qry_ret->aligned_seq.clear();
	qry_ret->insertions.clear();
	qry_ret->deletions.clear();
	qry_ret->substitutions.clear();
	qry_ret->match_pos.clear();
	//qry_ret->ret_seq=read->c_seq;
	int readLen=read->readLen;
	int OverlappingKmers=readLen-K_size+1;
	if(gap>=OverlappingKmers)
	{return;}
	int Read_arr_sz=readLen/32+1;
	int rem=readLen%32;
	if(rem==0)
	{Read_arr_sz--;}
	int tot_bits=Read_arr_sz*64;
	size_t ht_sz=ht->ht_sz;
	bool flip,found;
	size_t hash_idx;
	bool repeat_map=1,Smooth=0;
	
	
	
	uint64_t seq,f_seq,hv;
	bucket ** bktptr;
	char ref_char[500];
	char qry_char[500];

	for (int j=0;j<OverlappingKmers;j++)
	{

		get_sub_arr(read->read_bits,read->readLen,j,K_size,&(seq));	
		f_seq=get_rev_comp_seq(seq,K_size);
		flip=0;
		if(seq>f_seq)
		{
			uint64_t t=seq;
			seq=f_seq;
			f_seq=t;
			flip=1;
		}

//		if(key_map.count(seq))
	//	{map_match++;}

		hv=MurmurHash64A(&seq,sizeof(seq),0);

		hash_idx=(size_t) (hv%ht_sz);

		bktptr= &(ht->store_pos[hash_idx]);

		found=look_up_in_a_list(seq,&bktptr);
		if(found)
		{
			qry_ret->matching_kmers++;
			//flip=flip^((*bktptr)->kmer_info.flip);
			if(qry_ret->match_kmers_only)
			{
				for (int k=0;k<K_size;++k)
				{
					qry_ret->ret_seq[j+k]=toupper(qry_ret->ret_seq[j+k]);
				}
			
			}
		}
		else
		{continue;}
		
		if(1)
		{
			if((*bktptr)->kmer_info.repeat==0)
			{
				flip=flip^((*bktptr)->kmer_info.flip);
				if(flip==0)
				{
					//if((*bktptr)->kmer_info.cod>=j)
					//qry_ret->match_cnt[(*bktptr)->kmer_info.cod-j]++;
					
					qry_ret->match_cnt[(*bktptr)->kmer_info.cod]++;
					qry_ret->match_pos.push_back((*bktptr)->kmer_info.cod);
					
					//cout<<(*bktptr)->kmer_info.cod<<endl;
				}
				else
				{
					qry_ret->match_cnt[-(int64_t)((*bktptr)->kmer_info.cod)]++;
					qry_ret->match_pos.push_back(-(int64_t)((*bktptr)->kmer_info.cod));
					//qry_ret->match_cnt[-(int64_t)((*bktptr)->kmer_info.cod+K_size+j)]++;
					//cout<<-(int64_t)(*bktptr)->kmer_info.cod<<endl;
				}
			}
			else
			{
				if(repeat_map)
				{
					int map_idx=(*bktptr)->kmer_info.cod;
					if(ref->repeat_maps[map_idx].size()>=0xf)
					{
						continue;
					}
					else
					{
						for(map<uint64_t,bool>::iterator match_itr=(ref->repeat_maps[map_idx]).begin();match_itr!=ref->repeat_maps[map_idx].end();++match_itr)
						{
							bool flip2=flip;
							flip2=flip2^match_itr->second;
							if(flip2==0)
							{
								//cout<<int64_t(match_itr->first)-j<<endl;
								if(int64_t(match_itr->first)-j>=0)
								{qry_ret->match_cnt[int64_t(match_itr->first)-j]++;}
								else
								{qry_ret->match_cnt[0]++;}
							}
							else
							{
							
								//cout<<-(int64_t)(match_itr->first+K_size+j)<<endl;
								qry_ret->match_cnt[-(int64_t)(match_itr->first+K_size+j)]++;
							}
						
						}
					}
				}
			}
		}
	}

	if(!qry_ret->match_kmers_only)
	{
		if(qry_ret->match_cnt.size()==1)
		{
			//
			struct read_t read2;
			int64_t beg_pos=qry_ret->match_cnt.begin()->first;
			int band_width=2;
			if(beg_pos>=0)
			{
				int ref_len=read->readLen;
				if(beg_pos>=band_width&&(beg_pos+2*band_width<ref->ref_seq.readLen))
				{
					get_sub_arr(ref->ref_seq.read_bits,ref->ref_seq.readLen,beg_pos-band_width,read->readLen+2*band_width,read2.read_bits);
				
				}
				else
				{
					get_sub_arr(ref->ref_seq.read_bits,ref->ref_seq.readLen,beg_pos,read->readLen+2*band_width,read2.read_bits);
				}
				ref_len+=2*band_width;

				int Ref_arr_sz=ref_len/32+1;
				int rem=ref_len%32;
				if(rem==0)
				{Ref_arr_sz--;}
				bitsarr2str(read2.read_bits,ref_len,ref_char,Ref_arr_sz);
				qry_ret->matched_seq.push_back(ref_char);
				qry_ret->match_pos.push_back(beg_pos);
			}
			else
			{
				qry_ret->match_pos.push_back(beg_pos);
				beg_pos=-beg_pos-read->readLen;
			
				int ref_len=read->readLen;
				if(beg_pos>=band_width&&(beg_pos+2*band_width<ref->ref_seq.readLen))
				{
					get_sub_arr(ref->ref_seq.read_bits,ref->ref_seq.readLen,beg_pos-band_width,read->readLen+2*band_width,read2.read_bits);
				
				}
				else
				{
					get_sub_arr(ref->ref_seq.read_bits,ref->ref_seq.readLen,beg_pos,read->readLen,read2.read_bits);
				}
				ref_len+=2*band_width;
				int Ref_arr_sz=ref_len/32+1;
				int rem=ref_len%32;
				if(rem==0)
				{Ref_arr_sz--;}
				get_rev_comp_seq_arr(read2.read_bits,ref_len,Ref_arr_sz);
				bitsarr2str(read2.read_bits,ref_len,ref_char,Ref_arr_sz);
				qry_ret->matched_seq.push_back(ref_char);
			}
			if(qry_ret->align)
			{
				bitsarr2str(read->read_bits,read->readLen,qry_char,Read_arr_sz);
				//int match=100,mismatch=-100,gap_cost=-250,band_width=7;
				int match=100,mismatch=-100,gap_cost=-10,band_width=50;
				string A_aln,B_aln;
				struct aln_t aln_t; 
				//int b[300][300],score[300][300];
				int score=0; 
				//char qry_char2[300]="TCTTGCTTTTT";
				//char ref_char2[300]="TATTGCTTTTGAT";
			//	char qry_char2[300]="ATAGAAGTGGTGACAGAGGGCATCCTTGTCTTGTGCCAGTTTTCCACAATTTGCTTCCAGCTTTTGCTCATTCAGCATGATGTTAGCTTTGAGT";
			//	char ref_char2[300]="ATAGGAGTGGTGAGAGTGAGCATCTTTGTCATGTTTTTGTTCTCAAAGGGAATGCTTCCAGCTTTTGCTTGCTCACTATAATATTGGCTGTGGAT";
				//char qry_char2[300]="ATAGAAGTGGTGACAGAGGGCATCCTTGTCTTGTGCCAGTTTTCCACAATT      TGCTTCCAGCTTTTGCTCATTCAGCATGATGTTAGCTTTGAGT";
				//char ref_char2[300]="ATAGGAGTGGTGAGAGTGAGCATCTTTGTCATGTTTTTGTTCTCAAAGGGAA     TGCTTCCAGCTTTTGCTTGCTCACTATAATATTGGCTGTGGAT";
				score=GlobalAlign(&aln_t,qry_char,ref_char,match,mismatch, gap_cost,band_width);
				printAlign(&aln_t,qry_char,ref_char,A_aln,B_aln);
			//	cout<<A_aln<<endl;
			//	cout<<B_aln<<endl;


				//score=GlobalAlign(&aln_t,qry_char,ref_char,match,mismatch, gap_cost,band_width);
			
				printAlign(&aln_t,qry_char,ref_char,A_aln,B_aln);
				qry_ret->aligned_seq.push_back(A_aln);
				qry_ret->aligned_seq.push_back(B_aln);
				qry_ret->insertions.push_back(aln_t.insertions);
				qry_ret->deletions.push_back(aln_t.deletions);
				qry_ret->substitutions.push_back(aln_t.substitutions);
			}
	//		string ref_substr=
			//qry_ret->matched_seq.push_back();
		}
		else
		{
			if(qry_ret->match_cnt.size()==0)
			{return;}

			if(Smooth)
			{
				map<int64_t,int> max_match=qry_ret->match_cnt;
				map<int64_t,int>::iterator match_itr,p_match_itr,n_match_itr;
				p_match_itr=qry_ret->match_cnt.begin();
				for(match_itr=qry_ret->match_cnt.begin();match_itr!=qry_ret->match_cnt.end();)
				{
					n_match_itr=match_itr;
					n_match_itr++;
					if(n_match_itr==qry_ret->match_cnt.end())
					{break;}
					if(n_match_itr->first==match_itr->first+1)
					{
						max_match[match_itr->first]++;
					}
					if(p_match_itr->first==match_itr->first-1)
					{
						max_match[match_itr->first]++;
					}

					match_itr=n_match_itr;
					p_match_itr=match_itr;
				}
			}
			else
			{
				int64_t max_match_pos=0,max_cov=0;
				map<int64_t,int>::iterator match_itr;
				for(match_itr=qry_ret->match_cnt.begin();match_itr!=qry_ret->match_cnt.end();++match_itr)
				{
					if(match_itr->second>max_cov)
					{
						max_match_pos=match_itr->first;
						max_cov=match_itr->second;
					}
				}
				struct read_t read2;
				int64_t beg_pos=max_match_pos;
				int band_width=2;
				if(beg_pos>=0)
				{
				//	get_sub_arr(ref->ref_seq.read_bits,ref->ref_seq.readLen,beg_pos,read->readLen,read2.read_bits);
				//	bitsarr2str(read2.read_bits,read->readLen,ref_char,Read_arr_sz);

					int ref_len=read->readLen;
					if(beg_pos>=band_width&&(beg_pos+2*band_width<ref->ref_seq.readLen))
					{
					
						get_sub_arr(ref->ref_seq.read_bits,ref->ref_seq.readLen,beg_pos-band_width,read->readLen+2*band_width,read2.read_bits);
						//get_sub_arr(,)
					
					}
					else
					{
					
						get_sub_arr(ref->ref_seq.read_bits,ref->ref_seq.readLen,beg_pos,read->readLen+2*band_width,read2.read_bits);
					
					}
					ref_len+=2*band_width;
					int Ref_arr_sz=ref_len/32+1;
					int rem=ref_len%32;
					if(rem==0)
					{Ref_arr_sz--;}
				
					bitsarr2str(read2.read_bits,ref_len,ref_char,Ref_arr_sz);

					qry_ret->matched_seq.push_back(ref_char);
					qry_ret->match_pos.push_back(beg_pos);

				}
				else
				{
					qry_ret->match_pos.push_back(beg_pos);
					beg_pos=-beg_pos-read->readLen;


					int ref_len=read->readLen;
					if(beg_pos>=band_width&&(beg_pos+2*band_width<ref->ref_seq.readLen))
					{
						get_sub_arr(ref->ref_seq.read_bits,ref->ref_seq.readLen,beg_pos-band_width,read->readLen+2*band_width,read2.read_bits);
					
					}
					else
					{
						get_sub_arr(ref->ref_seq.read_bits,ref->ref_seq.readLen,beg_pos,read->readLen,read2.read_bits);
					}
					ref_len+=2*band_width;
					int Ref_arr_sz=ref_len/32+1;
					int rem=ref_len%32;
					if(rem==0)
					{Ref_arr_sz--;}

					get_rev_comp_seq_arr(read2.read_bits,ref_len,Ref_arr_sz);
					bitsarr2str(read2.read_bits,ref_len,ref_char,Ref_arr_sz);


					//get_sub_arr(ref->ref_seq.read_bits,ref->ref_seq.readLen,beg_pos,read->readLen,read2.read_bits);
					//get_rev_comp_seq_arr(read2.read_bits,read->readLen,Read_arr_sz);
					//bitsarr2str(read2.read_bits,read->readLen,ref_char,Read_arr_sz);
					qry_ret->matched_seq.push_back(ref_char);
				}
				if(qry_ret->align)
				{
					bitsarr2str(read->read_bits,read->readLen,qry_char,Read_arr_sz);
					//int match=100,mismatch=-100,gap_cost=-250,band_width=7;
					int match=100,mismatch=-100,gap_cost=-10,band_width=50;
					
					string A_aln,B_aln;
					struct aln_t aln_t; 
					int score=0; 
					score=GlobalAlign(&aln_t,qry_char,ref_char,match,mismatch, gap_cost,band_width);
					printAlign(&aln_t,qry_char,ref_char,A_aln,B_aln);
					qry_ret->aligned_seq.push_back(A_aln);
					qry_ret->aligned_seq.push_back(B_aln);
					qry_ret->insertions.push_back(aln_t.insertions);
					qry_ret->deletions.push_back(aln_t.deletions);
					qry_ret->substitutions.push_back(aln_t.substitutions);
				}

			}
		}

	}
	
}




void Sparse_Kmer_Ref_Graph_Query_LR(struct ref_t *ref, struct read_t *read,struct hashtable *ht,int K_size,int gap, struct qry_ret *qry_ret)
{
	qry_ret->matched_seq.clear();
	qry_ret->match_cnt.clear();
	qry_ret->matched_seq.clear();
	qry_ret->aligned_seq.clear();
	qry_ret->insertions.clear();
	qry_ret->deletions.clear();
	qry_ret->substitutions.clear();
	qry_ret->match_pos.clear();
	qry_ret->match_read_cod.clear();
	qry_ret->match_read.clear();
	qry_ret->ordered_read_pos.clear();
	qry_ret->match_stats.clear();
	qry_ret->match_ref_read_first.clear();
	qry_ret->match_ref_read_last.clear();

	int readLen=read->readLen;
	int OverlappingKmers=readLen-K_size+1;
	if(gap>=OverlappingKmers)
	{return;}
	int Read_arr_sz=readLen/32+1;
	int rem=readLen%32;
	if(rem==0)
	{Read_arr_sz--;}
	int tot_bits=Read_arr_sz*64;
	size_t ht_sz=ht->ht_sz;
	bool flip,found;
	size_t hash_idx;
	bool repeat_map=1,Smooth=0;
	
	uint64_t seq,f_seq,hv;
	bucket ** bktptr;
	char ref_char[500];
	char qry_char[500];
	map<int,vector<int> > read_position;
	for (int j=0;j<OverlappingKmers;j++)
	{

		get_sub_arr(read->read_bits,read->readLen,j,K_size,&(seq));	
		f_seq=get_rev_comp_seq(seq,K_size);
		flip=0;
		if(seq>f_seq)
		{
			uint64_t t=seq;
			seq=f_seq;
			f_seq=t;
			flip=1;
		}

		hv=MurmurHash64A(&seq,sizeof(seq),0);

		hash_idx=(size_t) (hv%ht_sz);

		bktptr= &(ht->store_pos[hash_idx]);

		found=look_up_in_a_list(seq,&bktptr);
		if(found)
		{
			qry_ret->matching_kmers++;
			if(qry_ret->match_kmers_only)
			{
				for (int k=0;k<K_size;++k)
				{
					qry_ret->ret_seq[j+k]=toupper(qry_ret->ret_seq[j+k]);
				}
			}
		}
		else
		{continue;}
		

		if((*bktptr)->kmer_info.repeat==0)
		{
			flip=flip^((*bktptr)->kmer_info.flip);
			vector <int64_t> read_pos_vt;
			int temp_read_idx=(*bktptr)->kmer_info.read_idx;

			read_pos_vt.push_back((*bktptr)->kmer_info.read_idx);//change here
			
			int read_pos=(*bktptr)->kmer_info.cod;
			int qry_len=ref->ref_read_len_vt[(*bktptr)->kmer_info.read_idx];
			if(flip==0)
			{
				qry_ret->match_read_cod[(*bktptr)->kmer_info.read_idx].push_back(read_pos);
				qry_ret->match_ref_read_beg[(*bktptr)->kmer_info.read_idx].push_back(j-read_pos);
				qry_ret->match_ref_read_end[(*bktptr)->kmer_info.read_idx].push_back(j+qry_len-read_pos);
				read_pos_vt.push_back((int32_t)(*bktptr)->kmer_info.cod);
				read_position[(*bktptr)->kmer_info.read_idx].push_back(read_pos);
			}
			else
			{
				qry_ret->match_read_cod[(*bktptr)->kmer_info.read_idx].push_back(-read_pos);
				qry_ret->match_ref_read_beg[-(int)(*bktptr)->kmer_info.read_idx].push_back(j+K_size-(qry_len-read_pos));
				qry_ret->match_ref_read_end[-(int)(*bktptr)->kmer_info.read_idx].push_back(j+read_pos+K_size);
				read_pos_vt.push_back(-(int32_t)(*bktptr)->kmer_info.cod);
				read_position[-(int)(*bktptr)->kmer_info.read_idx].push_back(read_pos);
			}

			int tmp_read_idx=(int)((*bktptr)->kmer_info.read_idx);
			if(flip==1)
			{
				tmp_read_idx=-tmp_read_idx;
			}

			if(qry_ret->match_ref_read_first.count(tmp_read_idx)==0)
			{
				qry_ret->match_ref_read_first[tmp_read_idx]=10000000;
			}
			if(qry_ret->match_ref_read_last.count(tmp_read_idx)==0)
			{
				qry_ret->match_ref_read_last[tmp_read_idx]=-10000000;
			}

			if(qry_ret->match_ref_read_first[tmp_read_idx]>j)
			{
				qry_ret->match_ref_read_first[tmp_read_idx]=j;
			}
			if(qry_ret->match_ref_read_last[tmp_read_idx]<j)
			{
				qry_ret->match_ref_read_last[tmp_read_idx]=j;
			}
			qry_ret->ordered_read_pos.push_back(read_pos_vt);
			
		}
		else
		{
			if(repeat_map)
			{
				

				int map_idx=(*bktptr)->kmer_info.cod;
				if(ref->repeat_maps[map_idx].size()>=0xff)
				{
					continue;
				}
				else
				{
					vector <int64_t> read_pos_vt;
					
					for(map<uint64_t,bool>::iterator match_itr=(ref->repeat_maps[map_idx]).begin();match_itr!=ref->repeat_maps[map_idx].end();++match_itr)
					{
						
						bool flip2=flip;
						flip2=flip2^match_itr->second;
						uint64_t cod_read=(match_itr->first);
						uint32_t read_idx= (uint32_t) cod_read;
						int tmp_read_idx=(int) read_idx;
						int qry_len=ref->ref_read_len_vt[read_idx];
						int read_pos;
						int32_t cod;
						if(flip2==0)
						{						
							cod= (int32_t) (cod_read>>32);
							read_pos=cod;
							qry_ret->match_ref_read_beg[tmp_read_idx].push_back(j-read_pos);
							qry_ret->match_ref_read_end[tmp_read_idx].push_back(j+qry_len-read_pos);
							read_position[tmp_read_idx].push_back(read_pos);
						}
						else
						{				
							tmp_read_idx=-tmp_read_idx;
							cod= -((int32_t) (cod_read>>32));
							read_pos=-cod;
							qry_ret->match_ref_read_beg[tmp_read_idx].push_back(j+K_size-(qry_len-read_pos));
							qry_ret->match_ref_read_end[tmp_read_idx].push_back(j+read_pos+K_size);
							read_position[tmp_read_idx].push_back(read_pos);
						}

						if(qry_ret->match_ref_read_first.count(tmp_read_idx)==0)
						{
							qry_ret->match_ref_read_first[tmp_read_idx]=1000000;
						}
						if(qry_ret->match_ref_read_last.count(tmp_read_idx)==0)
						{
							qry_ret->match_ref_read_last[tmp_read_idx]=-1000000;
						}

						if(qry_ret->match_ref_read_first[tmp_read_idx]>j)
						{
							qry_ret->match_ref_read_first[tmp_read_idx]=j;
						}
						if(qry_ret->match_ref_read_last[tmp_read_idx]<j)
						{
							qry_ret->match_ref_read_last[tmp_read_idx]=j;
						}

						qry_ret->match_read_cod[tmp_read_idx].push_back(cod);
						qry_ret->match_read.push_back(tmp_read_idx);
						qry_ret->match_pos.push_back(cod);
						read_pos_vt.push_back(tmp_read_idx);
						read_pos_vt.push_back(cod);
						
						
					}
					qry_ret->ordered_read_pos.push_back(read_pos_vt);
				}
				
			}
			
		}
	}

	for(map<int,vector<int> >::iterator tmp_itr=read_position.begin();tmp_itr!=read_position.end();++tmp_itr)
	{
		vector<int>::iterator vt_it1,tmp_it2;
		int min_pos=10000000,max_pos=-10000000;
		int match_bases=0,current_position=-100;
		int consecutive_matches=0,max_consecutive_matches=0;
		for(vt_it1=tmp_itr->second.begin();vt_it1!=tmp_itr->second.end();++vt_it1)
		{
			if(min_pos>*vt_it1)
			{min_pos=*vt_it1;}
			if(max_pos<*vt_it1)
			{max_pos=*vt_it1;}
			if(abs(*vt_it1-current_position)!=1)
			{
				match_bases+=K_size;
				if(consecutive_matches>max_consecutive_matches)
				{max_consecutive_matches=consecutive_matches;}
				consecutive_matches=0;

				
			}
			else
			{
				match_bases++;
				consecutive_matches++;
			}
			current_position=*vt_it1;

		}
		/*
		if(abs(max_pos-min_pos)<10&&match_bases>100)
		{
			for(vt_it1=tmp_itr->second.begin();vt_it1!=tmp_itr->second.end();++vt_it1)
			{
				cout<<*vt_it1<<endl;

			}
			cout<<endl;
		}
		*/
		
		qry_ret->match_stats[tmp_itr->first].push_back(min_pos);
		qry_ret->match_stats[tmp_itr->first].push_back(max_pos);
		qry_ret->match_stats[tmp_itr->first].push_back(match_bases);
		qry_ret->match_stats[tmp_itr->first].push_back(tmp_itr->second.size());
		qry_ret->match_stats[tmp_itr->first].push_back(max_consecutive_matches);
		//if(qry_ret->match_ref_read_first.count(abs(tmp_itr->first))==0)
		//{cout<<"error"<<endl;}
		//qry_ret->match_stats[tmp_itr->first].push_back(qry_ret->match_ref_read_first[abs(tmp_itr->first)]);
		//qry_ret->match_stats[tmp_itr->first].push_back(qry_ret->match_ref_read_last[abs(tmp_itr->first)]);
		
		qry_ret->match_stats[tmp_itr->first].push_back(qry_ret->match_ref_read_first[(tmp_itr->first)]);
		qry_ret->match_stats[tmp_itr->first].push_back(qry_ret->match_ref_read_last[(tmp_itr->first)]);
		
		
	}
	
	
}






/*

void Sparse_Kmer_Ref_Graph_Query_SSA(struct ref_t *ref, struct read_t *read,BF_info * BF_info,int K_size,int gap, struct qry_ret *qry_ret)
{
	qry_ret->matched_seq.clear();
	qry_ret->match_cnt.clear();
	qry_ret->matched_seq.clear();
	qry_ret->aligned_seq.clear();
	qry_ret->insertions.clear();
	qry_ret->deletions.clear();
	qry_ret->substitutions.clear();
	int readLen=read->readLen;
	int OverlappingKmers=readLen-K_size+1;
	if(gap>=OverlappingKmers)
	{return;}
	int Read_arr_sz=readLen/32+1;
	int rem=readLen%32;
	if(rem==0)
	{Read_arr_sz--;}
	int tot_bits=Read_arr_sz*64;
	size_t ht_sz=ht->ht_sz;
	bool flip,found;
	size_t hash_idx;
	bool repeat_map=1,Smooth=0;
	
	uint64_t seq,f_seq,hv;
	bucket ** bktptr;
	char ref_char[500];
	char qry_char[500];

	for (int j=0;j<OverlappingKmers;j++)
	{

		get_sub_arr(read->read_bits,read->readLen,j,K_size,&(seq));	
		f_seq=get_rev_comp_seq(seq,K_size);
		flip=0;
		if(seq>f_seq)
		{
			seq=f_seq;
			flip=1;
		}

		bf_inserted=1;

		for(int b=0;b<BF_info->d;++b)
		{
			size_t BF_hv=(size_t)(MurmurHash64A(&seq_c,sizeof(seq_c),b)%(BF_info->m));

			if((BF_info->BF_HT[BF_hv >> 3] & (1U << (BF_hv & 0x7U)))==0)
			{
				bf_inserted=0;
			}

		}
	
		if(bf_inserted)
		{
			first_found=j;
		}



		if((*bktptr)->kmer_info.repeat==0)
		{
			flip=flip^((*bktptr)->kmer_info.flip);
			if(flip==0)
			{
				if((*bktptr)->kmer_info.cod>=j)
				qry_ret->match_cnt[(*bktptr)->kmer_info.cod-j]++;
			}
			else
			{
				qry_ret->match_cnt[-(int64_t)((*bktptr)->kmer_info.cod+K_size+j)]++;
			}
		}

		
		else
		{
			if(repeat_map)
			{
				int map_idx=(*bktptr)->kmer_info.cod;
				if(ref->repeat_maps[map_idx].size()>=0xf)
				{
					continue;
				}
				else
				{
					for(map<uint64_t,bool>::iterator match_itr=(ref->repeat_maps[map_idx]).begin();match_itr!=ref->repeat_maps[map_idx].end();++match_itr)
					{
						bool flip2=flip;
						flip2=flip2^match_itr->second;
						if(flip2==0)
						{
							//cout<<int64_t(match_itr->first)-j<<endl;
							qry_ret->match_cnt[int64_t(match_itr->first)-j]++;
						}
						else
						{
							
							//cout<<-(int64_t)(match_itr->first+K_size+j)<<endl;
							qry_ret->match_cnt[-(int64_t)(match_itr->first+K_size+j)]++;
						}
						
					}
				}
			}
		}

	}
	if(qry_ret->match_cnt.size()==1)
	{
		//
		struct read_t read2;
		int64_t beg_pos=qry_ret->match_cnt.begin()->first;
		int band_width=2;
		if(beg_pos>=0)
		{
			int ref_len=read->readLen;
			if(beg_pos>=band_width&&(beg_pos+2*band_width<ref->ref_seq.readLen))
			{
				get_sub_arr(ref->ref_seq.read_bits,ref->ref_seq.readLen,beg_pos-band_width,read->readLen+2*band_width,read2.read_bits);
				
			}
			else
			{
				get_sub_arr(ref->ref_seq.read_bits,ref->ref_seq.readLen,beg_pos,read->readLen+2*band_width,read2.read_bits);
			}
			ref_len+=2*band_width;

			//must modify
			bitsarr2str(read2.read_bits,ref_len,ref_char,Read_arr_sz);
			qry_ret->matched_seq.push_back(ref_char);
			qry_ret->match_pos.push_back(beg_pos);
		}
		else
		{
			qry_ret->match_pos.push_back(beg_pos);
			beg_pos=-beg_pos-read->readLen;
			
			int ref_len=read->readLen;
			if(beg_pos>=band_width&&(beg_pos+2*band_width<ref->ref_seq.readLen))
			{
				get_sub_arr(ref->ref_seq.read_bits,ref->ref_seq.readLen,beg_pos-band_width,read->readLen+2*band_width,read2.read_bits);
				
			}
			else
			{
				get_sub_arr(ref->ref_seq.read_bits,ref->ref_seq.readLen,beg_pos,read->readLen,read2.read_bits);
			}
			ref_len+=2*band_width;
			//must modify
			get_rev_comp_seq_arr(read2.read_bits,ref_len,Read_arr_sz);
			bitsarr2str(read2.read_bits,ref_len,ref_char,Read_arr_sz);
			qry_ret->matched_seq.push_back(ref_char);
		}
		if(qry_ret->align)
		{
			bitsarr2str(read->read_bits,read->readLen,qry_char,Read_arr_sz);
			int match=100,mismatch=-100,gap_cost=-200,band_width=7;
			string A_aln,B_aln;
			struct aln_t aln_t; 
			//int b[300][300],score[300][300];
			int score=0; 
			//char qry_char2[300]="AAAACCCGTTTTT";
			//char ref_char2[300]="AAACCCCTTTTA";
		//	char qry_char2[300]="ATAGAAGTGGTGACAGAGGGCATCCTTGTCTTGTGCCAGTTTTCCACAATTTGCTTCCAGCTTTTGCTCATTCAGCATGATGTTAGCTTTGAGT";
		//	char ref_char2[300]="ATAGGAGTGGTGAGAGTGAGCATCTTTGTCATGTTTTTGTTCTCAAAGGGAATGCTTCCAGCTTTTGCTTGCTCACTATAATATTGGCTGTGGAT";
			//char qry_char2[300]="ATAGAAGTGGTGACAGAGGGCATCCTTGTCTTGTGCCAGTTTTCCACAATT      TGCTTCCAGCTTTTGCTCATTCAGCATGATGTTAGCTTTGAGT";
			//char ref_char2[300]="ATAGGAGTGGTGAGAGTGAGCATCTTTGTCATGTTTTTGTTCTCAAAGGGAA     TGCTTCCAGCTTTTGCTTGCTCACTATAATATTGGCTGTGGAT";
			score=GlobalAlign(&aln_t,qry_char,ref_char,match,mismatch, gap_cost,band_width);
			printAlign(&aln_t,qry_char,ref_char,A_aln,B_aln);
		//	cout<<A_aln<<endl;
		//	cout<<B_aln<<endl;


			//score=GlobalAlign(&aln_t,qry_char,ref_char,match,mismatch, gap_cost,band_width);
			
			printAlign(&aln_t,qry_char,ref_char,A_aln,B_aln);
			qry_ret->aligned_seq.push_back(A_aln);
			qry_ret->aligned_seq.push_back(B_aln);
			qry_ret->insertions.push_back(aln_t.insertions);
			qry_ret->deletions.push_back(aln_t.deletions);
			qry_ret->substitutions.push_back(aln_t.substitutions);
		}
//		string ref_substr=
		//qry_ret->matched_seq.push_back();
	}
	else
	{
		if(qry_ret->match_cnt.size()==0)
		{return;}

		if(Smooth)
		{
			map<int64_t,int> max_match=qry_ret->match_cnt;
			map<int64_t,int>::iterator match_itr,p_match_itr,n_match_itr;
			p_match_itr=qry_ret->match_cnt.begin();
			for(match_itr=qry_ret->match_cnt.begin();match_itr!=qry_ret->match_cnt.end();)
			{
				n_match_itr=match_itr;
				n_match_itr++;
				if(n_match_itr==qry_ret->match_cnt.end())
				{break;}
				if(n_match_itr->first==match_itr->first+1)
				{
					max_match[match_itr->first]++;
				}
				if(p_match_itr->first==match_itr->first-1)
				{
					max_match[match_itr->first]++;
				}

				match_itr=n_match_itr;
				p_match_itr=match_itr;
			}
		}
		else
		{
			int64_t max_match_pos=0,max_cov=0;
			map<int64_t,int>::iterator match_itr;
			for(match_itr=qry_ret->match_cnt.begin();match_itr!=qry_ret->match_cnt.end();++match_itr)
			{
				if(match_itr->second>max_cov)
				{
					max_match_pos=match_itr->first;
					max_cov=match_itr->second;
				}
			}
			struct read_t read2;
			int64_t beg_pos=max_match_pos;
			int band_width=2;
			if(beg_pos>=0)
			{
			//	get_sub_arr(ref->ref_seq.read_bits,ref->ref_seq.readLen,beg_pos,read->readLen,read2.read_bits);
			//	bitsarr2str(read2.read_bits,read->readLen,ref_char,Read_arr_sz);

				int ref_len=read->readLen;
				if(beg_pos>=band_width&&(beg_pos+2*band_width<ref->ref_seq.readLen))
				{
					get_sub_arr(ref->ref_seq.read_bits,ref->ref_seq.readLen,beg_pos-band_width,read->readLen+2*band_width,read2.read_bits);
					
				}
				else
				{
					get_sub_arr(ref->ref_seq.read_bits,ref->ref_seq.readLen,beg_pos,read->readLen+2*band_width,read2.read_bits);
				}
				ref_len+=2*band_width;
				//must modify
				bitsarr2str(read2.read_bits,ref_len,ref_char,Read_arr_sz);

				qry_ret->matched_seq.push_back(ref_char);
				qry_ret->match_pos.push_back(beg_pos);

			}
			else
			{
				qry_ret->match_pos.push_back(beg_pos);
				beg_pos=-beg_pos-read->readLen;


				int ref_len=read->readLen;
				if(beg_pos>=band_width&&(beg_pos+2*band_width<ref->ref_seq.readLen))
				{
					get_sub_arr(ref->ref_seq.read_bits,ref->ref_seq.readLen,beg_pos-band_width,read->readLen+2*band_width,read2.read_bits);
					
				}
				else
				{
					get_sub_arr(ref->ref_seq.read_bits,ref->ref_seq.readLen,beg_pos,read->readLen+2*band_width,read2.read_bits);
				}
				ref_len+=2*band_width;
				//must modify
				get_rev_comp_seq_arr(read2.read_bits,ref_len,Read_arr_sz);
				bitsarr2str(read2.read_bits,ref_len,ref_char,Read_arr_sz);


				//get_sub_arr(ref->ref_seq.read_bits,ref->ref_seq.readLen,beg_pos,read->readLen,read2.read_bits);
				//get_rev_comp_seq_arr(read2.read_bits,read->readLen,Read_arr_sz);
				//bitsarr2str(read2.read_bits,read->readLen,ref_char,Read_arr_sz);
				qry_ret->matched_seq.push_back(ref_char);
			}
			if(qry_ret->align)
			{
				bitsarr2str(read->read_bits,read->readLen,qry_char,Read_arr_sz);
				int match=100,mismatch=-100,gap_cost=-200,band_width=7;
				string A_aln,B_aln;
				struct aln_t aln_t; 
				int score=0; 
				score=GlobalAlign(&aln_t,qry_char,ref_char,match,mismatch, gap_cost,band_width);
				printAlign(&aln_t,qry_char,ref_char,A_aln,B_aln);
				qry_ret->aligned_seq.push_back(A_aln);
				qry_ret->aligned_seq.push_back(B_aln);
				qry_ret->insertions.push_back(aln_t.insertions);
				qry_ret->deletions.push_back(aln_t.deletions);
				qry_ret->substitutions.push_back(aln_t.substitutions);
			}

		}
	}

	
	
}

*/




#endif
