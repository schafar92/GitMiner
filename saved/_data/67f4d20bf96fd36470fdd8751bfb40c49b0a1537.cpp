#ifndef TAG_NODE_H
#define TAG_NODE_H
#include <iostream>
#include <sstream>
#include <chrono>
#include <numeric>
#include <limits>
#include <list>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <omp.h>


using namespace std;
//namespace HierarchicalClustering
//{
	class TagNode
	{
	public:
		//vector<int> tagIds;
		list<int> tagIds;
		TagNode* left;
		TagNode* right;
		float children_dist;
		//depth height
		//static shared_ptr<redis::client> conn;
		static unordered_map<int,string> tagDic;
		static unordered_map<int,unordered_map<int,float>> tagDist;

		TagNode():tagIds(),left(NULL),right(NULL),children_dist(0.0)
		{
		}

		TagNode(int tagId)
		{
			this->tagIds.push_back(tagId);
			this->left = NULL;
			this->right = NULL;
			this->children_dist = 0.0;
		}

		TagNode(list<int>& tagIds, TagNode* left,TagNode* right,float dist=0.0)
		{
			this->tagIds = tagIds;
			this->left = left;
			this->right = right;
			this->children_dist = dist;
		}

		TagNode(TagNode* left,TagNode* right,float dist=0.0)
		{
			this->tagIds.insert(this->tagIds.end(),left->tagIds.begin(),left->tagIds.end());
			this->tagIds.insert(this->tagIds.end(),right->tagIds.begin(),right->tagIds.end());
			
			this->left = left;
			this->right = right;
			this->children_dist = dist;
		}

		float distance(TagNode* other)
		{
			float sum = 0.0;
			for(auto i:this->tagIds)
			{
				for(auto j:other->tagIds)
				{
					sum += getTagDistance(i,j);
				}
			}
			return sum/(float)(this->tagIds.size()*other->tagIds.size());
		}

		int depth()
		{
			if(this->isLeaf())
			{
				return 0;
			}
			return max(this->left->depth(),this->right->depth()) + 1;
		}


		bool isLeaf()
		{
			return (this->left == NULL && this->right == NULL);
		}

		string toString()
		{
			stringstream ss;
			for(auto i:this->tagIds)
			{
				ss<<tagDic[i]<<" ";
			}
			return ss.str();
		}



		TagNode* merge(TagNode* other,float dist)
		{
			list<int> newVec;
			newVec.insert(newVec.end(),this->tagIds.begin(),this->tagIds.end());
			newVec.insert(newVec.end(),other->tagIds.begin(),other->tagIds.end());
			TagNode* newNode = new TagNode(newVec,this,other,dist);
			return newNode;

		}

		void query(int tagId, float threshold)
		{

		}
		
		float getTagDistance(int i,int j)
		{
			if(i<=0 || j<=0)
			{
				return 0.0;
			}
			else
			{
                                if(i>j)swap(i,j);
                        
				if(tagDist.find(i)!=tagDist.end())
				{
					return tagDist[i][j];
				}
				else
				{
					return 0.0;//return 1.0;
				}
			}
		}

	};

	unordered_map<int,string> TagNode::tagDic = unordered_map<int,string>();
	unordered_map<int,unordered_map<int,float>> TagNode::tagDist = unordered_map<int,unordered_map<int,float>>();

	TagNode* merge(TagNode* left,TagNode* right,float dist)
	{
		cout<<"max_distance:"<<dist<<endl;
		TagNode* newNode = new TagNode();
		list<int>& newVec = newNode->tagIds;
		newVec.insert(newVec.end(),left->tagIds.begin(),left->tagIds.end());
		newVec.insert(newVec.end(),right->tagIds.begin(),right->tagIds.end());
		newNode->left = left;
		newNode->right = right;
		newNode->children_dist = dist;
		return newNode;			
	}

	void loadData(vector<TagNode*>& nodes,string infile,string simfile)
	{
        const int max_w = 50;
        char tag[max_w];
        int tagId;
        FILE* f = fopen(infile.c_str(),"rb");
        if (f==NULL)
        {
                std::cout<<"errror reading dict file"<<std::endl;
                exit(1);
        }
        stringstream ss;
        while (fscanf(f,"%d\t%s ",&tagId,&tag )!=EOF)
        {       
            TagNode::tagDic[tagId]=string(tag);
            TagNode* np = new TagNode(tagId);
            nodes.push_back(np);
            //std::cout<<tagId<<" "<<tag<<std::endl;
        }
        fclose(f);
        tagId = 1205;//1205    做自己
        ss.clear();
        ss<<tagId;
        // string key = ss.str();
        // string s = TagNode::conn->get(key);

        std::cout<<"test retrive result:"<<endl;
        std::cout<<tagId<<" "<<TagNode::tagDic[tagId]<<std::endl;

        f = fopen(simfile.c_str(),"rb");
        if (f==NULL)
        {
                std::cout<<"errror reading dict file"<<std::endl;
                exit(1);
        }
        int tag1,tag2;
        float sim;
        auto start = std::chrono::system_clock::now();
		while (fscanf(f,"%d\t%d\t%f",&tag1,&tag2,&sim)!=EOF)
		{
	        if(tag1==tag2)continue;
	        if(tag1<tag2)
	        	TagNode::tagDist[tag1][tag2]=sim;
	        else
	        	TagNode::tagDist[tag2][tag1]=sim;

		}
		auto t = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now()-start).count();
		std::cout<<"loading tag similarity cost time: "<<t<<"ms"<<std::endl;
		fclose(f);

		int i = 2,j=4588;//2       1450    167
		cout<<i<<" "<<j<<" cosine similarity "<<TagNode::tagDist[i][j]<<endl;
    }

    TagNode* clustering(vector<TagNode*> Nodes)
    {
    	
    	while(Nodes.size()>1)
    	{
    		cout<<"clustering ["<<Nodes.size()<<"]"<<endl;
    		auto start = std::chrono::system_clock::now();
    		//float min_distance = numeric_limits<float>::max();
    		float max_distance = numeric_limits<float>::lowest();
    		int i=-1,j=-1;
    		int m,n;
			#pragma omp parallel
			{
			    float max_distance_private = max_distance;
			    int mi_private, mj_private;
			    #pragma omp for
			    for(int m=0;m<Nodes.size()-1;m++)
			    {
			        for(int n = m+1;n<Nodes.size();n++)
			        {
			            float distance = Nodes[m]->distance(Nodes[n]);
			            if(distance>max_distance_private)
			            {
			            	max_distance_private=distance;
			                mi_private=m;
			                mj_private=n;
			            }
			        }
			    }
			    #pragma omp critical 
			    {
			        if(max_distance_private>max_distance)
			        {
			            max_distance = max_distance_private;
			            i = mi_private;
			            j = mj_private;
			        }
			    }
			}

    		TagNode* node1 = Nodes[i];
    		TagNode* node2 = Nodes[j];
    		//delete j first because j>i
    		Nodes.erase(Nodes.begin()+j);
    		Nodes.erase(Nodes.begin()+i);
    		cout<<"merging node "<<i<<" and node "<<j<<":"<<endl;
    		// cout<<node1->toString()<<endl;
    		// cout<<node2->toString()<<endl;

    		TagNode* newNode = merge(node1,node2,max_distance);
    		Nodes.push_back(newNode);
			auto t = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now()-start).count();
			std::cout<<"this round of clustering cost: "<<t<<"ms"<<std::endl;
    	}
    	cout<<"clustering finished"<<endl;
    	return Nodes[0];
    }

    void clustering(vector<TagNode*> Nodes,float threshold,vector<TagNode*>& finalClusters)
    {
    	
    	while(Nodes.size()>1)
    	{
    		cout<<"clustering ["<<Nodes.size()<<"]"<<endl;
    		auto start = std::chrono::system_clock::now();
    		//float min_distance = numeric_limits<float>::max();
    		float max_distance = numeric_limits<float>::lowest();
    		int i=-1,j=-1;
    		int m,n;
			#pragma omp parallel
			{
			    float max_distance_private = max_distance;
			    int mi_private, mj_private;
			    #pragma omp for
			    for(int m=0;m<Nodes.size()-1;m++)
			    {
			        for(int n = m+1;n<Nodes.size();n++)
			        {
			            float distance = Nodes[m]->distance(Nodes[n]);
			            if(distance>max_distance_private)
			            {
			            	max_distance_private=distance;
			                mi_private=m;
			                mj_private=n;
			            }
			        }
			    }
			    #pragma omp critical 
			    {
			        if(max_distance_private>max_distance)
			        {
			            max_distance = max_distance_private;
			            i = mi_private;
			            j = mj_private;
			        }
			    }
			}

    		TagNode* node1 = Nodes[i];
    		TagNode* node2 = Nodes[j];
    		//delete j first because j>i
    		Nodes.erase(Nodes.begin()+j);
    		Nodes.erase(Nodes.begin()+i);
    		cout<<"merging node "<<i<<" and node "<<j<<":"<<endl;
    		// cout<<node1->toString()<<endl;
    		// cout<<node2->toString()<<endl;

    		TagNode* newNode = merge(node1,node2,max_distance);
    		Nodes.push_back(newNode);
			auto t = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now()-start).count();
			std::cout<<"this round of clustering cost: "<<t<<"ms"<<std::endl;
			if(max_distance<threshold)
				break;
    	}
    	cout<<"clustering finished"<<endl;
    	cout<<"filter cluster distance greater than threshold"<<endl;
    	//vector<TagNode*> finalClusters;
    	for(int m=0;m<Nodes.size();m++)
    	{
    		if(Nodes[m]->children_dist>threshold)
    		{
    			finalClusters.push_back(Nodes[m]);
    			cout<<Nodes[m]->toString()<<endl;
    		}
    	}
    	cout<<"final cluster size:"<<finalClusters.size()<<endl;

    	//return Nodes[0];
    }

    void levelOrder(TagNode* cur, int depth, vector<list<int>> &levelCluster,vector<float> &levelDist)
    {
    	if(cur==NULL || depth == 0)
    		return;
    	if(depth==1)
    	{
    		levelCluster.push_back(cur->tagIds);
    		levelDist.push_back(cur->children_dist);
    	}
    	depth-=1;
    	levelOrder(cur->left,depth,levelCluster,levelDist);
    	levelOrder(cur->right,depth,levelCluster,levelDist);

    }

    void cutTree(TagNode* root, float threshold)
    {
    	int depth = root->depth();
    	cout<<"tag hierarchy depth is "<<depth<<endl;
    	//vector<vector<list<int>>> levelClusters(depth,vector<list<int>>());
    	vector<vector<list<int>>> levelClusters;
    	int i;
       	for(i = 1;i<=depth;i++)
    	{
    		vector<list<int>> levelCluster;
    		vector<float> levelDist;
    		levelOrder(root,i,levelCluster,levelDist);
    		levelClusters.push_back(levelCluster);
    		bool cut = true;
    		//float ld = accumulate(levelDist.begin(),levelDist.end(),0.0);
    		//ld = ld/levelDist.size();
    		cout<<"level "<<i<<" and level size:"<<levelCluster.size()<<endl;
    		cout<<"level cluster distance:"<<endl;
    		for(auto ld:levelDist)
    		{
    			cout<<ld<<"\t";
    		}
    		cout<<endl;
    		for(auto j:levelDist)
    		{
    			if(j<threshold)
    			{
    				cut = false;
    				break;
    			}
    		}
    		if(cut)
    		{
    			break;
    		}
    	}
    	cout<<"cut level "<<i<<" found and cluster size:"<<levelClusters[i].size()<<endl;
    	for(auto levelCluster:levelClusters[i])
    	{
    		for(auto tagId:levelCluster)
    		{
    			cout<<TagNode::tagDic[tagId]<<"\t";
    		}
    		cout<<endl;
    	}
    }

    string neighbors(TagNode* node, int queryId)
    {
		stringstream ss;
		for(auto i:node->tagIds)
		{
			if(i!=queryId)
				ss<<TagNode::tagDic[i]<<" ";
		}
		return ss.str();
    }

    void query(vector<TagNode*>& clusters,int queryId)
    {
    	TagNode* target=NULL;
    	for(auto cluster:clusters)
    	{
    		if(find(cluster->tagIds.begin(),cluster->tagIds.end(),queryId)!=cluster->tagIds.end())
    		//if(cluster->tagIds.find(queryId)!=cluster->tagIds.end())
    		{
    			//cout<<cluster->toString()<<endl;
    			target = cluster;
    		}
    	}
    	if(target)
    	{
    		TagNode* node = target;
    		int level = 0;
    		while(node)
    		{
    			if(node->tagIds.size()<=1)
    				break;
    			//cout<<"level:"<<level++<<"\t"<<node->toString()<<endl;
    			cout<<"level:"<<level++<<"\t"<<neighbors(node,queryId)<<endl;
    			if(node->left)
    			{
    				list<int>& cur = node->left->tagIds;
    				if(find(cur.begin(),cur.end(),queryId)!=cur.end())
    				//if(node->left->tagIds.find(queryId)!=node->left->tagIds.end())
    				{
    					node = node->left;
    					continue;
    				}
    			}
    			if(node->right)
    			{
    				list<int>& cur = node->right->tagIds;
    				if(find(cur.begin(),cur.end(),queryId)!=cur.end())
    				//if(node->right->tagIds.find(queryId)!=node->right->tagIds.end())
    				{
    					node = node->right;
    					continue;
    				}
    			}
    			else
    				break;
    		}
    		
    	}
    	else
    	{
    		cout<<queryId<<" not found in the tag clusters"<<endl;
    	}
    }

	int main(int argc, char* argv[])
	{
        TagNode a(0);
        cout<<sizeof(a)<<endl;
		if(argc<2)
		{
			cout<<"usage:exe dic_file sim_file"<<endl;
		}
		string dicfile(argv[1]);
		string simfile(argv[2]);
		vector<TagNode*> tagNodes,finalClusters;
        loadData(tagNodes,dicfile,simfile);
		//TagNode* root = clustering(tagNodes);
		float threshold = 0.40;
		clustering(tagNodes,threshold,finalClusters);
		int queryId = 10074;
		query(finalClusters,queryId);
		//cut the tag hierarchy to tree
		//float threshold = 0.60;
		//cutTree(root,threshold);
        return 0;
	}
//}
#endif
