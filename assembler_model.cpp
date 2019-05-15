#include <bits/stdc++.h>
#include <fstream>
#include <vector>
#include <string>

#include "tensorflow/cc/ops/const_op.h"
#include "tensorflow/cc/ops/image_ops.h"
#include "tensorflow/cc/ops/standard_ops.h"
#include "tensorflow/core/framework/graph.pb.h"
#include "tensorflow/core/framework/tensor.h"
#include "tensorflow/core/graph/default_device.h"
#include "tensorflow/core/graph/graph_def_builder.h"
#include "tensorflow/core/lib/core/errors.h"
#include "tensorflow/core/lib/core/stringpiece.h"
#include "tensorflow/core/lib/core/threadpool.h"
#include "tensorflow/core/lib/io/path.h"
#include "tensorflow/core/lib/strings/stringprintf.h"
#include "tensorflow/core/platform/init_main.h"
#include "tensorflow/core/platform/logging.h"
#include "tensorflow/core/platform/types.h"
#include "tensorflow/core/public/session.h"
#include "tensorflow/core/util/command_line_flags.h"

#define pb push_back
#define eb emplace_back
#define ll long long
using namespace std;

struct kmer{
	int utgno;
	vector<int>rid;
};

struct contig{
	int ctind;
	string s;
	set<int>sid;
	int nxt;
	int prv;
	int hdr;
}cntg;

int k,abd,iter2,idx,ress;
float mx_pairedEnd;

vector<contig>ctgs,ctgsn;
vector<string>edg;
vector<int>vtmp;
set<int>st,tip,hdrst;

unordered_map<string,kmer>hmap,hmap2;
unordered_map<string,kmer>::iterator it;
set<int>::iterator itt;

ofstream opf("unitigs.fasta");
ofstream opf2("unitig_info.csv");
ofstream opf3("contig.fasta");
ofstream opf4("branch_info.csv");
ofstream opf5("kmer_count.txt");

string cmap[4]={"A","T","G","C"};
string orig,tmp;

bool cmp(contig a,contig b){
	if(a.sid.size()==b.sid.size())
		return a.s.length()>b.s.length();
	return a.sid.size()>b.sid.size();
}

void generateUnitigs(unordered_map<string,kmer>::iterator it){
	orig=(*it).first;
	cntg.s=orig;
	cntg.sid.clear();
	idx=ctgs.size();
	tmp=orig;
	while(true){
		hmap[tmp].utgno=idx;
		edg.clear();
		for(int i=0;i<hmap[tmp].rid.size();i++){
			cntg.sid.insert(hmap[tmp].rid[i]);
		}
		hmap[tmp].rid.clear();
		for(int i=0;i<4;i++){
			if(hmap.find(tmp.substr(1,k)+cmap[i])!=hmap.end()){ 
				edg.pb(cmap[i]);
			}
		}
		if(edg.size()!=1)
			break;
		else{
			tmp=tmp.substr(1,k)+edg[0];
			if(hmap[tmp].rid.size()<1)
				break;
			for(int i=0;i<4;i++){
				if(hmap.find(cmap[i]+tmp.substr(0,k-1))!=hmap.end()){ 
					edg.pb(cmap[i]);
				}
			}
			if(edg.size()==2){
				cntg.s+=tmp[k-1];
			}
			else{
				break;
			}
		}
	}
	tmp=orig;
	while(true){
		hmap[tmp].utgno=idx;
		edg.clear();
		for(int i=0;i<hmap[tmp].rid.size();i++){
			cntg.sid.insert(hmap[tmp].rid[i]);
		}
		hmap[tmp].rid.clear();
		for(int i=0;i<4;i++){
			if(hmap.find(cmap[i]+tmp.substr(0,k-1))!=hmap.end()){ 
				edg.pb(cmap[i]);
			}
		}
		if(edg.size()!=1){
			break;
		}
		else{
			tmp=edg[0]+tmp.substr(0,k-1);
			if(hmap[tmp].rid.size()<1)
				break;
			for(int i=0;i<4;i++){
				if(hmap.find(tmp.substr(1,k)+cmap[i])!=hmap.end()){ 
					edg.pb(cmap[i]);
				}
			}
			if(edg.size()==2){
				cntg.s=edg[0]+cntg.s;
			}
			else{
				break;
			}
		}
	}
	iter2++;
	if(iter2%5000==0)
		cout<<iter2<<"\n";
	opf<<">contig "<<iter2<<"\n";
	opf<<cntg.s<<"\n";
	cntg.nxt=0;
	cntg.prv=0;
	cntg.ctind=ctgs.size();
	cntg.hdr=cntg.ctind;
	hdrst.insert(cntg.ctind);
	ctgs.eb(cntg);
	return ;
}

void make_csv(){
	opf2<<"unitig,readid_of_unitig,right_end_A,readid,right_end_T,readid,right_end_G,readid,right_end_C,readid,left_end_A,readid,left_end_T,readid,left_end_G,readid,left_end_C,readid\n";
	int idx;
	for(int i=1;i<ctgs.size();i++){
		vtmp.clear();
		int cntr=0,cntl=0;
		string strr=ctgs[i].s.substr(ctgs[i].s.length()-k+1,k-1);
		string strl=ctgs[i].s.substr(0,k-1);
		for(int j=0;j<4;j++){
			if((hmap.find(strr+cmap[j])!=hmap.end()) && (hmap[strr+cmap[j]].utgno!=i) && (hmap[strr+cmap[j]].utgno>0)){
				vtmp.eb(hmap[strr+cmap[j]].utgno);
				cntr++;
			}
			else{
				vtmp.eb(-1);
			}
		}
		for(int j=0;j<4;j++){
			if((hmap.find(cmap[j]+strl)!=hmap.end()) && (hmap[cmap[j]+strl].utgno!=i) && (hmap[cmap[j]+strl].utgno>0)){
				vtmp.eb(hmap[cmap[j]+strl].utgno);
				cntl++;
			}
			else{
				vtmp.eb(-1);
			}
		}
		if(cntr>1 || cntl>1){
			opf2<<ctgs[i].s<<",";
			for(itt=ctgs[i].sid.begin();itt!=ctgs[i].sid.end();itt++)
				opf2<<(*itt)<<"#";
			opf2<<",";
			for(int j=0;j<8;j++){
				if(vtmp[j]==-1)
					opf2<<"$,$,";
				else{
					opf2<<ctgs[vtmp[j]].s<<",";
					for(itt=ctgs[vtmp[j]].sid.begin();itt!=ctgs[vtmp[j]].sid.end();itt++)
						opf2<<(*itt)<<"#";
					opf2<<",";
				}
			}
			opf2<<"\n";
		}		
	}
}


void make_csv2(){
	set<int>tmpp;
	int tpcnt=0;
	int vecc[5];
	bool mrkk[100000];
	fill(mrkk,mrkk+100000,false);
	opf4<<"A,T,G,C\n";
	for(int i=1;i<ctgs.size();i++){
		string suff=ctgs[i].s.substr(ctgs[i].s.length()-k+1,k-1);
		fill(vecc,vecc+5,0);
		tpcnt=0;
		for(int j=0;j<4;j++){
			tmpp.clear();
			if(hmap.find(suff+cmap[j])!=hmap.end()){
				int indd=hmap[suff+cmap[j]].utgno;
				for(int x=0;x<k-1;x++){
					for(int y=0;y<hmap2[ctgs[indd].s.substr(x,k)].rid.size();y++){
						tmpp.insert(hmap2[ctgs[indd].s.substr(x,k)].rid[y]);
					}
				}	
				/*for(itt=ctgs[indd].sid.begin();itt!=ctgs[indd].sid.end();itt++){
					if((*itt)&1){
						if(ctgs[i].sid.find((*itt)+1)!=ctgs[i].sid.end())
							tpcnt++;
					}
					else{
						if(ctgs[i].sid.find((*itt)-1)!=ctgs[i].sid.end())
							tpcnt++;
					}
				}*/
				mrkk[indd]=true;	
			}
			vecc[j]=tmpp.size();
			if(tmpp.size())
				tpcnt++;
		}
		if(tpcnt>1){
			for(int j=0;j<3;j++)
				opf4<<vecc[j]<<",";
			opf4<<vecc[3]<<"\n";
		}
	}
	for(int i=1;i<ctgs.size();i++){
		if(mrkk[i]==false)continue;
		tmpp.clear();
		for(int x=0;x<k-1;x++){
			for(int y=0;y<hmap2[ctgs[i].s.substr(x,k)].rid.size();y++){
				tmpp.insert(hmap2[ctgs[i].s.substr(x,k)].rid[y]);
			}
		}
		opf5<<tmpp.size()<<"\n";
	}
}



void tips_removal(){
	int l=0,r=0;
	for(int i=1;i<ctgs.size();i++){
		string suff=ctgs[i].s.substr(ctgs[i].s.length()-k+1,k-1);
		string pre=ctgs[i].s.substr(0,k-1);
		r=0;l=0;
		for(int j=0;j<4;j++){
			if(hmap.find(suff+cmap[j])!=hmap.end())
				r++;
			if(hmap.find(cmap[j]+pre)!=hmap.end())
				l++;
		}
		if((l==0 || r==0) && ctgs[i].s.length()<2*k)
			tip.insert(i);
	}
	cout<<"no. of tips "<<tip.size()<<"\n";
}
	
//the machine learning model used to for branch resolution. It is not working properly.
/*int evaluate_model(string contig)
{

	Session* session;

	Status status = NewSession(SessionOptions(), &session);
	if (!status.ok()) {
	std::cout << status.ToString() << "\n";
	return 1;
	}
	std::string graph_definition = "./model/mlp.pb";

	GraphDef graph_def;
	SessionOptions opts;

	std::vector<Tensor> outputs; 
	TF_CHECK_OK(ReadBinaryProto(Env::Default(), graph_definition, &graph_def));


	graph::SetDefaultDevice("/gpu:0", &graph_def);
	opts.config.mutable_gpu_options()->set_per_process_gpu_memory_fraction(0.5);
	opts.config.mutable_gpu_options()->set_allow_growth(true);


	TF_CHECK_OK(NewSession(opts, &session));
	TF_CHECK_OK(session->Create(graph_def));


	TF_CHECK_OK(session->Run({}, {}, {"init_all_vars_op"}, nullptr));

	Tensor x(DT_FLOAT, TensorShape({1, 11})); // number of features (contig length, repeats etc...)
	Tensor y(DT_FLOAT, TensorShape({1, 6}));  // number of output labels

	// auto _XTensor = x.matrix<float>();
	// auto _YTensor = y.matrix<float>();

	std::vector<Tensor> outputs;
	string input_layer = "Mul";
	string output_layer = "softmax";
	if (!session->Run({{input_layer, x}},
	                 {output_layer}, {}, &outputs).ok()) {
	LOG(ERROR) << "Running model failed";
	return -1;
	}

	TF_CHECK_OK(session.Run({{x, {data_set.input(get_features (contig) )}}}, {layer_3}, &outputs));
	Eigen::Map<Eigen::VectorXf> pred(outputs[0].flat<float>().data(),
	                               outputs[0].NumElements());
	int maxIndex; float maxValue = pred.maxCoeff(&maxIndex);
	session->Close();


	return maxIndex;
  }



}*/



void branching(){
	int cntr=0,cntr2=0;
	int cntt[100000],cntt2[100000];
	fill(cntt,cntt+100000,0);
	fill(cntt2,cntt2+100000,0);
	set<int>tmpp;
	for(int i=1;i<ctgs.size();i++)
		ctgsn.eb(ctgs[i]);
	sort(ctgsn.begin(),ctgsn.end(),cmp);
	int mxx,mxind,iter3=0;
	string suff;
	int m_score, cum_score;
	for(int i=0;i<ctgsn.size();i++){
		int ind=ctgsn[i].ctind;
		if(ctgs[ind].prv!=0 || ind==0 || (tip.find(ind)!=tip.end()))continue;	
		while(ctgs[ind].nxt==0){
			mxx=0;mxind=0;
			suff=ctgs[ind].s.substr(ctgs[ind].s.length()-k+1,k-1);
			for(int j=0;j<4;j++){
				int tpcnt=0;
				if(hmap.find(suff+cmap[j])!=hmap.end()){
					if(ctgs[hmap[suff+cmap[j]].utgno].prv==0 && (tip.find(hmap[suff+cmap[j]].utgno)==tip.end())){
						int indd=hmap[suff+cmap[j]].utgno;
						
						for(itt=ctgs[ctgs[indd].hdr].sid.begin();itt!=ctgs[ctgs[indd].hdr].sid.end();itt++){
							if((*itt)&1){
								if(ctgs[ctgs[ind].hdr].sid.find((*itt)+1)!=ctgs[ctgs[ind].hdr].sid.end())
									tpcnt++;
							}
							else{
								if(ctgs[ctgs[ind].hdr].sid.find((*itt)-1)!=ctgs[ctgs[ind].hdr].sid.end())
									tpcnt++;
							}
						}

						// Use model when length of contig is greater than 500
						/*if ((ctgs[ctgs[ind].hdr].s.length() + ctgs[ctgs[indd].hdr].s.length()) > 500 )
						{
							m_score = evaluate_model(ctgs[ctgs[indd].hdr].s);
						}
						else
						{
							m_score =-1;
						}*/
						if(tpcnt>=0 && m_score!=-1)
						{
							if (tpcnt>mx_pairedEnd)
								mx_pairedEnd= tpcnt;
							tpcnt/= mx_pairedEnd
							m_score /= 5.0;
							cum_score = tpcnt + m_score;
						}
						/*if (m_score!= -1)
						{
							cum_score = m_score;
						}*/
						if(cum_score>mxx){	//check if this branch's cumulative score is the highest, if yes set this branch as mxind; 
							mxx=cum_score;
							mxind=indd;
						}
					}
				}
			}
			/*if(mxx==0){
				for(int j=0;j<4;j++){
					if(hmap.find(suff+cmap[j])!=hmap.end()){
						if(ctgs[hmap[suff+cmap[j]].utgno].prv==0 && (tip.find(hmap[suff+cmap[j]].utgno)==tip.end())){
							tmpp.clear();
							int indd=hmap[suff+cmap[j]].utgno;
							for(int x=0;x<k-1;x++){
								for(int y=0;y<hmap2[ctgs[indd].s.substr(x,k)].rid.size();y++){
									tmpp.insert(hmap2[ctgs[indd].s.substr(x,k)].rid[y]);
									cntr++;
								}
							}
							if(tmpp.size()>mxx ){
								mxx=tmpp.size();
								mxind=hmap[suff+cmap[j]].utgno;
							}
							//if(ctgs[hmap[suff+cmap[j]].utgno].s.length()>mxx){
								//mxx=ctgs[hmap[suff+cmap[j]].utgno].s.length();
								//mxind=hmap[suff+cmap[j]].utgno;
							//}
						}
					}
				}
			}*/
			if(mxind==0)break;
			ctgs[ind].nxt=mxind;
			ctgs[mxind].prv=ind;
			for(itt=ctgs[ctgs[mxind].hdr].sid.begin();itt!=ctgs[ctgs[mxind].hdr].sid.end();itt++)
				ctgs[ctgs[ind].hdr].sid.insert(*itt);
			ctgs[ctgs[ind].hdr].s+=ctgs[ctgs[mxind].hdr].s;
			hdrst.erase(ctgs[mxind].hdr);
			ctgs[mxind].hdr=ctgs[ind].hdr;
			ind=mxind;
		}
	}
	for(itt=hdrst.begin();itt!=hdrst.end();itt++){
		if((*itt)==0||(tip.find(*itt)!=tip.end()))continue;
		opf3<<">contig "<<++iter3<<"\n";
		opf3<<ctgs[*itt].s<<"\n";
	}
}


					
		

int main(int argc, char** argv){
	ios_base::sync_with_stdio(false);
	auto start=chrono::system_clock::now();
	int iter=0;
	k=stoi(argv[2]);
	abd=stoi(argv[3]);
	freopen(argv[1],"r",stdin);
	
	string read,str;
	while(getline(cin,read)){
		getline(cin,read);
		// preprocessing to replace missassemblies with 'C'
		for(int i=0;i<read.length();i++){
			if(read[i]!='A' && read[i]!='T' && read[i]!='G' && read[i]!='C')
				read[i]='C';
		}
		for(int i=0;i+k<=read.length();i++){
			str=read.substr(i,k);
			hmap[str].rid.emplace_back(iter);
				
		}
		getline(cin,read);
		getline(cin,read);
		if(iter%100000==0)
			cout<<iter<<"\n";
		iter++;
	}
	cntg.s="$";
	ctgs.eb(cntg);
	cout<<"Stage 1 done\n";
	vector<unordered_map<string,kmer>::iterator>tmpp;
	cout<<"ddd "<<hmap.size()<<endl;
	for(it=hmap.begin();it!=hmap.end();it++){
		hmap2[(*it).first]=(*it).second;
		if((*it).second.rid.size()<abd){
			tmpp.eb(it);
		}
	}
	for(int i=0;i<tmpp.size();i++)
		hmap.erase(tmpp[i]);
	tmpp.clear();
	cout<<"ddd "<<hmap.size()<<endl;
	for(it=hmap.begin();it!=hmap.end();it++){
		if((*it).second.rid.size())
			generateUnitigs(it);
	}
	cout<<"Stage 2 done\n";
	//make_csv();
	//tips_removal();
	branching();
	//make_csv2();
	cout<<"Stage 3 done\n";
	cout<<"ddd "<<hmap.size()<<endl;
	auto end=chrono::system_clock::now();
	chrono::duration<double>el=end-start;
	cout<<el.count();
	return 0;
}
