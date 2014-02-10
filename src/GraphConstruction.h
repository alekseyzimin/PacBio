#ifndef __GRAPH_CONSTRUCTION_H
#define __GRAPH_CONSTRUCTION_H


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


// initialize a hashtable
void Init_HT(struct hashtable* ht,size_t ht_sz)
{
	ht->ht_sz=ht_sz;
	ht->store_pos=(struct bucket**)calloc(ht_sz,sizeof(struct bucket*));

	for(size_t i=0;i<ht_sz;++i)
	{
		ht->store_pos[i]=NULL;
	}
}

void Init_HT2(struct hashtable2* ht,size_t ht_sz)
{
	ht->ht_sz=ht_sz;
	ht->store_pos=(struct bucket2**)calloc(ht_sz,sizeof(struct bucket2*));
	for(size_t i=0;i<ht_sz;++i)
	{
		ht->store_pos[i]=NULL;
	}
}

void Init_HT3(struct hashtable3* ht,size_t ht_sz)
{
	ht->ht_sz=ht_sz;
	ht->store_pos=(struct bucket3**)calloc(ht_sz,sizeof(struct bucket3*));
	for(size_t i=0;i<ht_sz;++i)
	{
		ht->store_pos[i]=NULL;
	}
}

void Init_HT4(struct hashtable4* ht,size_t ht_sz)
{
	ht->ht_sz=ht_sz;
	ht->store_pos=(struct bucket4**)calloc(ht_sz,sizeof(struct bucket4*));
	for(size_t i=0;i<ht_sz;++i)
	{
		ht->store_pos[i]=NULL;
	}
}






//look up for a k-mer in a hashtable, if exists: 1, otherwise: 0. for round 2
bool look_up_in_a_list(uint64_t seq,struct bucket *** ptr)
{
	bool found=0;
	//struct bucket ** bktptr;


	while((**ptr)!=NULL)
	{
		if((**ptr)->kmer_t.kmer==seq)
		{
		//	bktptr=ptr;
			break;
		}

		(*ptr)=&((**ptr)->nxt_bucket);
	}
	if((**ptr)==NULL)
	{
		found=0;
		//bktptr=ptr;
	}
	else
	{
		found=1;
	}
	return found;
}

bool look_up_in_a_list2(struct kmer_t2 *seq,struct bucket2 *** ptr)
{
	bool found=0;
	//struct bucket ** bktptr;

	while((**ptr)!=NULL)
	{
		if(memcmp(&((**ptr)->kmer_t2.kmer),&(seq->kmer),sizeof(uint64_t)*2)==0)
		{
			break;
		}

		(*ptr)=&((**ptr)->nxt_bucket);
	}
	if((**ptr)==NULL)
	{
		found=0;
	}
	else
	{
		found=1;
	}
	return found;
}

bool look_up_in_a_list3(struct kmer_t3 *seq,struct bucket3 *** ptr)
{
	bool found=0;
	//struct bucket ** bktptr;

	while((**ptr)!=NULL)
	{
		if(memcmp(&((**ptr)->kmer_t3.kmer),&(seq->kmer),sizeof(uint64_t)*3)==0)
		{
			break;
		}

		(*ptr)=&((**ptr)->nxt_bucket);
	}
	if((**ptr)==NULL)
	{
		found=0;
	}
	else
	{
		found=1;
	}
	return found;
}

bool look_up_in_a_list4(struct kmer_t4 *seq,struct bucket4 *** ptr)
{
	bool found=0;
	//struct bucket ** bktptr;

	while((**ptr)!=NULL)
	{
		if(memcmp(&((**ptr)->kmer_t4.kmer),&(seq->kmer),sizeof(uint64_t)*4)==0)
		{
			break;
		}

		(*ptr)=&((**ptr)->nxt_bucket);
	}
	if((**ptr)==NULL)
	{
		found=0;
	}
	else
	{
		found=1;
	}
	return found;
}


bool look_up_in_a_list_rm(uint64_t seq,struct bucket_rm *** ptr)
{
	bool found=0;
	//struct bucket ** bktptr;


	while((**ptr)!=NULL)
	{
		if((**ptr)->kmer_t.kmer==seq)
		{
		//	bktptr=ptr;
			break;
		}

		(*ptr)=&((**ptr)->nxt_bucket);
	}
	if((**ptr)==NULL)
	{
		found=0;
		//bktptr=ptr;
	}
	else
	{
		found=1;
	}
	return found;
}



bool look_up_in_a_list_rm2(struct kmer_t2 *seq,struct bucket_rm2 *** ptr)
{
	bool found=0;
	//struct bucket ** bktptr;

	while((**ptr)!=NULL)
	{
		if(memcmp(&((**ptr)->kmer_t2.kmer),&(seq->kmer),sizeof(uint64_t)*2)==0)
		{
			break;
		}

		(*ptr)=&((**ptr)->nxt_bucket);
	}
	if((**ptr)==NULL)
	{
		found=0;
	}
	else
	{
		found=1;
	}
	return found;
}

bool look_up_in_a_list_rm3(struct kmer_t3 *seq,struct bucket_rm3 *** ptr)
{
	bool found=0;
	//struct bucket ** bktptr;

	while((**ptr)!=NULL)
	{
		if(memcmp(&((**ptr)->kmer_t3.kmer),&(seq->kmer),sizeof(uint64_t)*3)==0)
		{
			break;
		}

		(*ptr)=&((**ptr)->nxt_bucket);
	}
	if((**ptr)==NULL)
	{
		found=0;
	}
	else
	{
		found=1;
	}
	return found;
}

bool look_up_in_a_list_rm4(struct kmer_t4 *seq,struct bucket_rm4 *** ptr)
{
	bool found=0;
	//struct bucket ** bktptr;

	while((**ptr)!=NULL)
	{
		if(memcmp(&((**ptr)->kmer_t4.kmer),&(seq->kmer),sizeof(uint64_t)*4)==0)
		{
			break;
		}

		(*ptr)=&((**ptr)->nxt_bucket);
	}
	if((**ptr)==NULL)
	{
		found=0;
	}
	else
	{
		found=1;
	}
	return found;
}



// compare 2 arrays.
int uint64_t_cmp(uint64_t* A,uint64_t* B,int Kmer_arr_sz)
{
	int flag=0;
	for(int jj=0;jj<Kmer_arr_sz;++jj)
	{
		if(A[jj]>B[jj])
		{
			flag=1;
			break;
		}
		if(A[jj]<B[jj])
		{
			flag=-1;
			break;
		}
		if(A[jj]==B[jj])
		{
			continue;
		}
	}
	return flag;


}


//look up for a k-mer in a hashtable, if exists: 1, otherwise: 0. for round 1
bool look_up_in_a_list_r1(uint64_t seq,struct bucket_r1 *** ptr)
{	
	bool found=0;

	while((**ptr)!=NULL)
	{
		if((**ptr)->kmer_t.kmer==seq)
		{
		//	bktptr=ptr;	
			break;
		}

		(*ptr)=&((**ptr)->nxt_bucket);
	}
	if((**ptr)==NULL)
	{
		found=0;
		//bktptr=ptr;
	}
	else
	{
		found=1;
	}
	return found;
}

bool look_up_in_a_list2_r1(struct kmer_t2 *seq,struct bucket2_r1 *** ptr)
{	
	bool found=0;
	
	while((**ptr)!=NULL)
	{
		if(memcmp(&((**ptr)->kmer_t2.kmer),&(seq->kmer),sizeof(uint64_t)*2)==0)
		{
			break;
		}

		(*ptr)=&((**ptr)->nxt_bucket);
	}
	if((**ptr)==NULL)
	{
		found=0;
	}
	else
	{
		found=1;
	}
	return found;
}

bool look_up_in_a_list3_r1(struct kmer_t3 *seq,struct bucket3_r1 *** ptr)
{	
	bool found=0;
	
	while((**ptr)!=NULL)
	{
		if(memcmp(&((**ptr)->kmer_t3.kmer),&(seq->kmer),sizeof(uint64_t)*3)==0)
		{
			break;
		}

		(*ptr)=&((**ptr)->nxt_bucket);
	}
	if((**ptr)==NULL)
	{
		found=0;
	}
	else
	{
		found=1;
	}
	return found;
}

bool look_up_in_a_list4_r1(struct kmer_t4 *seq,struct bucket4_r1 *** ptr)
{	
	bool found=0;
	
	while((**ptr)!=NULL)
	{
		if(memcmp(&((**ptr)->kmer_t4.kmer),&(seq->kmer),sizeof(uint64_t)*4)==0)
		{
			break;
		}

		(*ptr)=&((**ptr)->nxt_bucket);
	}
	if((**ptr)==NULL)
	{
		found=0;
	}
	else
	{
		found=1;
	}
	return found;
}




//reference graph



void Sparse_Kmer_Ref_Graph_Construction(struct ref_read_t *read,struct hashtable *ht,int64_t *bucket_count,int64_t *edge_cnt,int K_size,int gap,BF_info * BF_info ,int round,int ref_pos,struct ref_t *ref)//,map<uint64_t,int > &key_map)
{
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
	bool flip[200],found[200];
	size_t hash_idx[200];
	


//	map<uint64_t,int > key_map;
/*
	for(int j=0;j<OverlappingKmers;++j)
	{
		uint64_t test_bits=0,test_bits2=0;	
		get_sub_arr(read->read_bits,read->readLen,j,K_size,&test_bits);	
		test_bits2=get_rev_comp_seq(test_bits,K_size);
		if(test_bits>test_bits2)
		{
			test_bits=test_bits2;
		}
		key_map[test_bits]++;

	}
	cout<<"bucket_count1:"<< key_map.size()<<endl;

	*/

	memset(flip,0,sizeof(flip));
	memset(found,0,sizeof(found));
	uint64_t seq[200],f_seq[200],hv[200];
	bucket ** bktptr[200];
	bool found_t=0,found_c=0,found_n=0,flip_t=0,flip_c=0,flip_n=0;
	uint64_t seq_c,seq_n,f_seq_n,hv_n;
	size_t hash_idx_t;
	bucket ** bktptr_c,** bktptr_n;
	
	int first_found=-1;
	int pre_search_len=100;
	if(OverlappingKmers<pre_search_len)
	{pre_search_len=OverlappingKmers;}
	for (int j=0;j<pre_search_len;j++)
	{
		
		get_sub_arr(read->read_bits,read->readLen,j,K_size,&(seq[j]));	
		f_seq[j]=get_rev_comp_seq(seq[j],K_size);
		if(seq[j]>f_seq[j])
		{
			uint64_t t=seq[j];
			seq[j]=f_seq[j];
			f_seq[j]=t;
			flip[j]=1;
		}
	

		hv[j]=MurmurHash64A(&seq[j],sizeof(seq[j]),0);

		hash_idx[j]=(size_t) (hv[j]%ht_sz);

		bktptr[j]= &(ht->store_pos[hash_idx[j]]);

		if(round==1)
		{
			found[j]=look_up_in_a_list_r1(seq[j], ((bucket_r1 ***) &bktptr[j] ));
			if(found[j]==1)
			{
				if(first_found<0)
				{first_found=j;}
			}
		}
		else
		{
			found[j]=look_up_in_a_list(seq[j],&bktptr[j]);
		}

	}
	int g,h;

	g=0;
	for (int k=0;k<gap;++k)
	{
		if(round==1)
		{
			found[k]=look_up_in_a_list_r1(seq[k],(bucket_r1 ***)&bktptr[k]);
		}
		if(found[k]==1)
		{
			g=k;
			
			break;
		}

	}

	//check if there is a saved one in the BF table
	
	bool bf_inserted=0;
	if(round==1&&found[g]==0&&BF_info->Bloom)
	{
		
		for (int k=0;k<gap;++k)
		{
			bf_inserted=1;
			if(BF_info->Bloom) 
			{
				
				for(int b=0;b<BF_info->d;++b)
				{
					size_t BF_hv=(size_t)(MurmurHash64A(&seq[k],sizeof(seq[k]),b)%(BF_info->m));
	
					if((BF_info->BF_HT[BF_hv >> 3] & (1U << (BF_hv & 0x7U)))==0)
					{
						bf_inserted=0;
					}
					//BF_info->BF_HT[BF_hv >> 3] |= (1U << (BF_hv & 0x7U));
				}
			}
			if(bf_inserted)
			{
				g=k;break;
			}

		}
	}
	
	///////updated for better sparseness
	if(round==1&&first_found>0&&found[g]==0)
	{
		g=first_found%gap;

	}
	found_c=found[g];
	seq_c=seq[g];
	bktptr_c=bktptr[g];
	flip_c=flip[g];
	//look for the next sparse kmer

	for(int j=g;j<OverlappingKmers;)
	{
		
		

	
		h=gap;
		bool bf_inserted=1;
		for (int k=1;k<=gap;++k)
		{
			
			//bf_inserted=1;	
			if( (j+k) == OverlappingKmers )//overbound
			{
				h=k+1;
				break;
			}
			//update the found
			
			if(j+k<pre_search_len)
			{
				if(round==1)
				{
					
					if(BF_info->Bloom)
					{
						
						bf_inserted=1;
						for(int b=0;b<BF_info->d;++b)
						{
							size_t BF_hv=(size_t)(MurmurHash64A(&seq[j+k],sizeof(seq[j+k]),b)%(BF_info->m));

							if((BF_info->BF_HT[BF_hv >> 3] & (1U << (BF_hv & 0x7U)))==0)
							{
								bf_inserted=0;
								//BF_info->BF_HT[BF_hv >> 3] |= (1U << (BF_hv & 0x7U));
								
							}
					
						}
						
					}
					found_n=0;
					if((!BF_info->Bloom)||(BF_info->Bloom&&bf_inserted))
					{
						
						found_n=look_up_in_a_list_r1(seq[j+k],(bucket_r1 ***) &bktptr[j+k]);
					}
				}
				else
				{
					found_n=look_up_in_a_list(seq[j+k],&bktptr[j+k]);
				}
				bktptr_n=bktptr[j+k];
				flip_n=flip[j+k];
				seq_n=seq[j+k];
			}
			else
			{
				get_sub_arr(read->read_bits,read->readLen,j+k,K_size,&(seq_n));	
				f_seq_n=get_rev_comp_seq(seq_n,K_size);
				flip_n=0;
				if(seq_n>f_seq_n)
				{
					uint64_t t=seq_n;
					seq_n=f_seq_n;
					f_seq_n=t;
					flip_n=1;
				}

				hv_n=MurmurHash64A(&seq_n,sizeof(seq_n),0);

				hash_idx_t=(size_t) (hv_n%ht_sz);

				bktptr_n= &(ht->store_pos[hash_idx_t]);

				if(round==1)
				{
					if(BF_info->Bloom)
					{
						
						bf_inserted=1;
						for(int b=0;b<BF_info->d;++b)
						{
							size_t BF_hv=(size_t)(MurmurHash64A(&seq_n,sizeof(seq_n),b)%(BF_info->m));

							if((BF_info->BF_HT[BF_hv >> 3] & (1U << (BF_hv & 0x7U)))==0)
							{
								bf_inserted=0;
								//BF_info->BF_HT[BF_hv >> 3] |= (1U << (BF_hv & 0x7U));
								
							}
					
						}
						
					}
					found_n=0;
					if((!BF_info->Bloom)||(BF_info->Bloom&&(bf_inserted==1)))
					{
						
					found_n=look_up_in_a_list_r1(seq_n,(bucket_r1 ***) &bktptr_n);
					}
				}
				else
				{
					found_n=look_up_in_a_list(seq_n,&bktptr_n);
				}
			}
			

			if(found_n==1)
			{
				h=k;
				
				break;
			}

		}

		
		if(round==1)
		{
			if(BF_info->Bloom)
			{
				
				bf_inserted=1;
				for(int b=0;b<BF_info->d;++b)
				{
					size_t BF_hv=(size_t)(MurmurHash64A(&seq_c,sizeof(seq_c),b)%(BF_info->m));

					if((BF_info->BF_HT[BF_hv >> 3] & (1U << (BF_hv & 0x7U)))==0)
					{
						bf_inserted=0;
						BF_info->BF_HT[BF_hv >> 3] |= (1U << (BF_hv & 0x7U));
								
					}
					
				}
						
			}
			if((!BF_info->Bloom)||(BF_info->Bloom&&(bf_inserted==1)))
			{
				
				found_c=look_up_in_a_list_r1(seq_c,(bucket_r1 ***) &bktptr_c);
			}
		}
		if(round==2)
		{
			found_c=look_up_in_a_list(seq_c,&bktptr_c);
		}
		if(found_c==0)
		{
			if(round==1)
			{
				*(bktptr_c)=(struct bucket*)malloc(sizeof(struct bucket_r1));
				memset(*bktptr_c,0,sizeof(struct bucket_r1));
				//memset((ptr->kmer_t),0,sizeof(struct ptr->kmer_t));
				//ptr=ptr->nxt_bucket;
				((struct bucket_r1*) *(bktptr_c))->nxt_bucket=NULL;
				((struct bucket_r1*) *(bktptr_c))->kmer_t.kmer=seq_c;
				//(*(bktptr[j]))->kmer_info.cov1=0;
				((struct bucket_r1*) *(bktptr_c))->kmer_info.cov1++;
				//(*(bktptr[j]))->kmer_info.left=0;
				//(*(bktptr[j]))->kmer_info.right=0;
				//(*(bktptr[j]))->kmer_info.split_left=0;
				//(*(bktptr[j]))->kmer_info.split_right=0;
				//(*(bktptr[j]))->kmer_info.used=0;
				//(*(bktptr[j]))->kmer_info.cod=0;
				//(*(bktptr[j]))->kmer_info.contig_no=0;
				(*bucket_count)++;
			}
		}
		else
		{
			if(round==1)
			{
				if(((struct bucket_r1*) *(bktptr_c))->kmer_info.cov1<0xffff)
				{((struct bucket_r1*) *(bktptr_c))->kmer_info.cov1++;}
			}
			else
			{
				if((*(bktptr_c))->kmer_info.cov1<0xffff)
				{
					(*(bktptr_c))->kmer_info.cov1++;
					if((*(bktptr_c))->kmer_info.cov1==1)
					{
						(*(bktptr_c))->kmer_info.cod=j;//ref_pos+
						(*(bktptr_c))->kmer_info.read_idx=read->read_idx;
						(*(bktptr_c))->kmer_info.flip=flip_c;
					}
					else
					{

						if((*(bktptr_c))->kmer_info.cov1==2)
						{
							(*(bktptr_c))->kmer_info.repeat=1;
							uint64_t tmp_cod=(*(bktptr_c))->kmer_info.cod;
							bool tmp_flip=(*(bktptr_c))->kmer_info.flip;
							(*(bktptr_c))->kmer_info.cod=ref->repeat_cnt;
							(ref->repeat_cnt)++;
							map<uint64_t,bool> r_map;
							//cout<<"cod: "<<tmp_cod<<endl;
							tmp_cod<<=32;//higher 32 bits are set for cod
							tmp_cod|=(*(bktptr_c))->kmer_info.read_idx;//lower 32 bits are set for read_idx
							r_map[tmp_cod]=tmp_flip;
							ref->repeat_maps.push_back(r_map);

						}
						size_t map_idx=((*(bktptr_c))->kmer_info.cod);
						uint64_t tmp_cod=j;//ref_pos+
						tmp_cod<<=32;
						tmp_cod|=read->read_idx;
						ref->repeat_maps[map_idx][tmp_cod]=flip_c;

						
					
					}
					
				}
			}
		}


		found_c=found_n;
		if(seq_c==seq_n)
		{
			found_c=1;
		}
		bktptr_c=bktptr_n;
		seq_c=seq_n;
		flip_c=flip_n;
		j=j+h;

	}


}




//convert the bucket type from round 1 to round 2. The buckets in round 2 are more expensive

void SwitchBuckets(hashtable *ht,hashtable2 *ht2,int K_size)
{
	size_t ht_sz;
	if(K_size<=32)
	{
		ht_sz=ht->ht_sz;
		bucket_r1 *store_pos_o,*store_pos_t;
		bucket *store_pos_n;
		bucket **bktp2p;
		for(size_t i=0;i<ht_sz;++i)
		{
			bktp2p=&(ht->store_pos[i]);
			store_pos_o=(bucket_r1*) ht->store_pos[i];
			while(store_pos_o!=NULL)
			{
				store_pos_n=(bucket*) malloc(sizeof(struct bucket));
				memset(store_pos_n,0,sizeof(struct bucket));
				store_pos_n->kmer_t=store_pos_o->kmer_t;
				store_pos_n->kmer_info.cov1=store_pos_o->kmer_info.cov1;
				store_pos_n->kmer_info.left=NULL;
				store_pos_n->kmer_info.right=NULL;
				*bktp2p=store_pos_n;
				bktp2p=&(store_pos_n->nxt_bucket);
				store_pos_t=store_pos_o;
				store_pos_o=store_pos_o->nxt_bucket;
				free(store_pos_t);
			}
		}
	}

	if(K_size>32&&K_size<64)
	{
		ht_sz=ht2->ht_sz;
		bucket2_r1 *store_pos_o,*store_pos_t;
		bucket2 *store_pos_n;
		bucket2 **bktp2p;
		for(size_t i=0;i<ht_sz;++i)
		{
			bktp2p=&(ht2->store_pos[i]);
			store_pos_o=(bucket2_r1*) ht2->store_pos[i];

			/*
			
			int n_buckets=0;
			while(store_pos_o!=NULL)
			{
				n_buckets++;
				store_pos_o=store_pos_o->nxt_bucket;
			
			}
			store_pos_o=(bucket2_r1*) ht2->store_pos[i];
			store_pos_n=(bucket2*) malloc(sizeof(struct bucket2)*n_buckets);
			bucket2 * init_bucket=store_pos_n;
			memset(store_pos_n,0,sizeof(bucket2)*n_buckets);
			n_buckets=0;	
			while(store_pos_o!=NULL)
			{
				store_pos_n->kmer_t2=store_pos_o->kmer_t2;
				store_pos_n->kmer_info.cov1=store_pos_o->kmer_info.cov1;
				store_pos_n->kmer_info.left=NULL;
				store_pos_n->kmer_info.right=NULL;
				*bktp2p=store_pos_n;
				bktp2p=&(store_pos_n->nxt_bucket);
				store_pos_t=store_pos_o;
				store_pos_o=store_pos_o->nxt_bucket;
				free(store_pos_t);
				//store_pos_n+=sizeof(bucket2);
				n_buckets++;
				store_pos_n=&(init_bucket[n_buckets]);
				
			}
			*/


			
			while(store_pos_o!=NULL)
			{
				store_pos_n=(bucket2*) malloc(sizeof(struct bucket2));
				memset(store_pos_n,0,sizeof(bucket2));
				store_pos_n->kmer_t2=store_pos_o->kmer_t2;
				store_pos_n->kmer_info.cov1=store_pos_o->kmer_info.cov1;
				store_pos_n->kmer_info.left=NULL;
				store_pos_n->kmer_info.right=NULL;
				*bktp2p=store_pos_n;
				bktp2p=&(store_pos_n->nxt_bucket);
				store_pos_t=store_pos_o;
				store_pos_o=store_pos_o->nxt_bucket;
				free(store_pos_t);
			}
			
		}
	}

}












#endif
