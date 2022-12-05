#include<iostream>
#include<fstream>
#include<unordered_map>
#include<map>
#include<unordered_set>
#include<queue>
#include<set>
#include<limits.h>
#include<vector>
#include <sys/stat.h>
using namespace std;

#define BLOCKSIZE 2
#define K 10

// #define NODES "testNodes.txt"
// #define EDGES "testEdges.txt"

#define NODES "sample_nodes.txt"
#define EDGES "sample_edges.txt"

// #define NODES "nodes.txt"
// #define EDGES "edges.txt"

struct point{
	long double x, y;
	string block_x_y;
	point(long double x, long double y, string block_x_y){
		this->x = x;
		this->y = y;
		this->block_x_y = block_x_y;
	}
};

struct edge{
	int u, v;
	long double length;
	edge(int u, int v, long double length){
		this->u = u;
		this->v = v;
		this->length = length;
	}
};

unordered_map<int, struct point*> nodeid_cordinates;

void insert(string fileName, unordered_map<int, struct point*> &nodeid_cordinates, vector<int> &nodes, vector<int> &border_nodes, unordered_set<struct edge*> &inedge, unordered_set<struct edge*> &outedge ){

	// Write internal nodes to file
	int num = nodes.size();
	mkdir("data", 0777);
	ofstream ofp("data/"+fileName+".txt");
	ofp << "----------\n";
	int count = 0;
	int overflowNumber = 0;
	for(int i=0; i<num; i++){
		int id = nodes[i];
		long double x = nodeid_cordinates[id]->x;
		long double y = nodeid_cordinates[id]->y;
		ofp << setprecision (numeric_limits<double>::digits10 + 1) << id << "  " << "<" << x << ", " << y << ">\n";
		count++;
		if(count == BLOCKSIZE){
			overflowNumber++;
			count = 0;
			string overflowFile = fileName+"_"+to_string(overflowNumber)+".txt";
			ofp << "?? "+overflowFile << "\n";
			ofp.close();
			ofp.open("data/"+overflowFile);
			// ofp << "\n?? "+fileName+".txt" << "\n";
			ofp << "?? "+fileName+".txt" << "\n";
		}
	}

	// write internal edges into file
	ofp << "##\n";
	for(unordered_set<struct edge*>:: iterator i=inedge.begin(); i!=inedge.end(); i++){
		int u = (*i)->u;
		int v = (*i)->v;
		long double length = (*i)->length;
		ofp << setprecision (numeric_limits<double>::digits10 + 1) << u << " " << v << " <" << length << ">\n";
		count++;
		if(count == BLOCKSIZE){
			overflowNumber++;
			count = 0;
			string overflowFile = fileName+"_"+to_string(overflowNumber)+".txt";
			ofp << "?? "+overflowFile << "\n";
			ofp.close();
			ofp.open("data/"+overflowFile);
			ofp << "?? "+fileName+".txt" << "\n";
		}
	}

	// write border nodes to file
	ofp << "**\n";
	for(int i=0; i<border_nodes.size(); i++){
		int id = border_nodes[i];
		long double x = nodeid_cordinates[id]->x;
		long double y = nodeid_cordinates[id]->y;
		ofp << setprecision (numeric_limits<double>::digits10 + 1) << id << "  " << "<" << x << ", " << y << ">\n";
		count++;
		if(count == BLOCKSIZE){
			overflowNumber++;
			count = 0;
			string overflowFile = fileName+"_"+to_string(overflowNumber)+".txt";
			ofp << "?? "+overflowFile << "\n";
			ofp.close();
			ofp.open("data/"+overflowFile);
			ofp << "?? "+fileName+".txt" << "\n";
		}
	}

	// write border egdes to file
	ofp << "%%\n";
	int outEdgeCount = outedge.size();
	for(unordered_set<struct edge*>:: iterator i=outedge.begin(); i!=outedge.end(); i++){
		int u = (*i)->u;
		int v = (*i)->v;
		long double length = (*i)->length;
		ofp << setprecision (numeric_limits<double>::digits10 + 1) << u << " " << v << " <" << length << ">\n";
		count++;
		outEdgeCount--;
		if(count == BLOCKSIZE && outEdgeCount > 0){
			overflowNumber++;
			count = 0;
			string overflowFile = fileName+"_"+to_string(overflowNumber)+".txt";
			ofp << "?? "+overflowFile << "\n";
			ofp.close();
			ofp.open("data/"+overflowFile);
			ofp << "?? "+fileName+".txt" << "\n";
		}
	}
	ofp << "----------\n";
	ofp.close();
}

void createPartition(unordered_map<int, struct point*> &nodeid_cordinates, unordered_map<string, vector<int> > &block_node, long double x_min, long double y_min){
	for(unordered_map<int, struct point*>::iterator itr = nodeid_cordinates.begin(); itr != nodeid_cordinates.end(); itr++){
		int nodeid = itr->first;
		int x_rem = (itr->second->x - x_min) / K;
		int y_rem = (itr->second->y - y_min) / K;
		string s = to_string(x_rem) + "," + to_string(y_rem);
		itr->second->block_x_y = s;
		block_node[s].push_back(nodeid);
	}
}

void insertIntoDataBlocks(unordered_map<int, struct point*> &nodeid_cordinates, unordered_map<string, vector<int> > &block_node, unordered_map<string, vector<int> > &border_nodes, unordered_map<string, unordered_set<struct edge*> > &block_inedge, unordered_map<string, unordered_set<struct edge*> > &block_outedge){
	for(unordered_map<string, vector<int> >::iterator itr = block_node.begin(); itr != block_node.end(); itr++){
		string block_x_y = itr->first;
		vector<int> nodes = itr->second;
		vector<int> border_nodes_list = border_nodes[block_x_y];
		unordered_set<struct edge*> inedge = block_inedge[block_x_y];
		unordered_set<struct edge*> outedge = block_outedge[block_x_y];
		insert(block_x_y, nodeid_cordinates, nodes, border_nodes_list, inedge, outedge);
	}
}

int partitionHelperFunction(){
	unordered_map<string, vector<int> > block_node, border_nodes;
	unordered_map<string, unordered_set<struct edge*> > block_inedge, block_outedge;

	//Parse nodes file
	ifstream ifpn(NODES);
	if(!ifpn) return -1;
	int nodeid;
	long double x_min = INT_MAX, y_min = INT_MAX, x_max = INT_MIN, y_max = INT_MIN;
	while(ifpn >> nodeid){
		long double temp1, temp2;
		ifpn >> temp1 >> temp2;

		//map: nide_id -> (x, y, cell_id)
		struct point *p = new point(temp1, temp2, "");
		nodeid_cordinates[nodeid] = p;

		// find the x_min, y_min, x_max, y_max
		if(temp1 < x_min) x_min = temp1;
		else if(temp1 > x_max) x_max = temp1;

		if(temp2 < y_min) y_min = temp2;
		else if(temp2 > y_max) y_max = temp2;
	}
	ifpn.close();

	createPartition(nodeid_cordinates, block_node, x_min, y_min);

	//Parse edges file
	ifstream ifpe(EDGES);
	unordered_set<string> edgesDuplicates;
	if(!ifpe) return -1;
	int u, v;
	while(ifpe >> u){
		long double len;
		ifpe >> v;
		ifpe >> len;

		//to remove duplicate edges
		string temp1 = to_string(u)+"_"+to_string(v);
		string temp2 = to_string(v)+"_"+to_string(u);
		if(edgesDuplicates.find(temp1) != edgesDuplicates.end()) continue;
		edgesDuplicates.insert(temp1);
		edgesDuplicates.insert(temp2);

		struct edge* e = new edge(u, v, len);
		string block_u = nodeid_cordinates[u]->block_x_y;
		string block_v = nodeid_cordinates[v]->block_x_y;
		if(block_u == block_v){
			block_inedge[block_u].insert(e);
		} else {
			block_outedge[block_u].insert(e);
			block_outedge[block_v].insert(e);
			if(find(border_nodes[block_u].begin(), border_nodes[block_u].end(), v) == border_nodes[block_u].end()){
				border_nodes[block_u].push_back(v);
			}
			if(find(border_nodes[block_v].begin(), border_nodes[block_v].end(), v) == border_nodes[block_v].end()){
				border_nodes[block_v].push_back(u);
			}
		}
	}
	edgesDuplicates.clear();
	ifpe.close();

	//insert into file
	insertIntoDataBlocks(nodeid_cordinates, block_node, border_nodes, block_inedge, block_outedge);

	return 0;
}

int findNodeId(string record){
	string temp;
	for(int i=0; i<record.length() && record[i] != '\t' && record[i] != ' '; i++){
		temp += record[i];
	}
	try{
		return stoi(temp);
	} catch(...){
		cout << temp << "\n";
		return -1;
	}
}

void findNodeUV(string record, int *u, int *v, double *w){
	string temp1, temp2, temp3;
	int i=0;
	for(; i<record.length() && record[i] != '\t' && record[i] != ' '; i++){
		temp1 += record[i];
	}
	i++;
	for(; i<record.length() && record[i] != '\t' && record[i] != ' '; i++){
		temp2 += record[i];
	}
	i += 2;
	for(; i<record.length() && record[i] != '\t' && record[i] != '>'; i++){
		temp3 += record[i];
	}
	try{
		*u = stoi(temp1);
		*v = stoi(temp2);
		*w = stod(temp3);
	} catch(...){
		cout << temp1 << " ";
		cout << temp2 << " ";
		cout << temp3 << " ";
		cout << "\n";
	}
}

map<int, long double> d;
map<int, bool> visited;
map<int, int> parent;
map<int, vector<pair<int, long double> > > adj;
set<string> filesRead;
int blockCount = 0;

void readFromFile(string blockName, int overflowNumber, bool inode, bool iedge, bool bnode, bool bedge){
	ifstream ifp;
	string fileName = overflowNumber == 0 ? blockName+".txt": blockName+"_"+to_string(overflowNumber)+".txt";
	fileName = "data/"+fileName;
	ifp.open(fileName);
	if(!ifp){
		return;
	}
	filesRead.insert(fileName);
	blockCount++;
	// Read internal nodes
	string temp;
	getline (ifp, temp); // read file line
	if(inode){
		while (!ifp.eof() && getline (ifp, temp)) {
			if(temp[0] == '?'){
				overflowNumber++;
				return readFromFile(blockName, overflowNumber, inode, iedge, bnode, bedge);
			}
			if(temp == "##"){inode = false; break;}
			int nodeID = findNodeId(temp);
			if(nodeID == -1) cout << fileName << "  1\n";
			if(d.find(nodeID) == d.end()) d[nodeID] = INT_MAX;
			if(visited.find(nodeID) == visited.end()) visited[nodeID] = false;
			if(parent.find(nodeID) == parent.end()) parent[nodeID] = -1;
		}
	}
	// Read internal edges
	if(iedge){
		while (!ifp.eof() && getline (ifp, temp)) {
			if(temp[0] == '?'){
				overflowNumber++;
				return readFromFile(blockName, overflowNumber, inode, iedge, bnode, bedge);
			}
			if(temp == "**"){iedge = false; break;}
			int u, v;
			double w;
			findNodeUV(temp, &u, &v, &w);
			pair<int, long double> vw, uw;
	    vw.first = v; vw.second = w; uw.first = u; uw.second = w;
	    adj[u].push_back(vw);
	    adj[v].push_back(uw);
		}
	}
	// Read external nodes
	if(bnode){
		while (!ifp.eof() && getline (ifp, temp)) {
			if(temp[0] == '?'){
				overflowNumber++;
				return readFromFile(blockName, overflowNumber, inode, iedge, bnode, bedge);
			}
			if(temp == "%%"){bnode = false; break;}
			int nodeID = findNodeId(temp);
			if(nodeID == -1) cout << fileName << "  2";
			if(d.find(nodeID) == d.end()) d[nodeID] = INT_MAX;
		}
	}
	// Read external edges
	if(bedge){
		while (!ifp.eof() && getline (ifp, temp)) {
			if(temp[0] == '?'){
				overflowNumber++;
				return readFromFile(blockName, overflowNumber, inode, iedge, bnode, bedge);
			}
			if(temp == "----------"){bedge = false; return;}
			int u, v;
			double w;
			findNodeUV(temp, &u, &v, &w);
			pair<int, long double> vw, uw;
	    vw.first = v; vw.second = w; uw.first = u; uw.second = w;
	    adj[u].push_back(vw);
	    adj[v].push_back(uw);
		}
	}
	overflowNumber++;
	readFromFile(blockName, overflowNumber, inode, iedge, bnode, bedge);
}

void readFromFileHelperFunction(string blockName){
	string fileName = blockName + ".txt";
	if(filesRead.find(fileName) == filesRead.end()){
		readFromFile(blockName, 0, true, true, true, true);
	}
}

class Compare{
public:
  int operator() (edge* v1, edge* v2){
    return v1->length > v2->length;
  }
};


void printPath(int dst){
	vector<int> result;
	while(parent[dst] != -1){
		result.push_back(dst);
		dst = parent[dst];
	}
	result.push_back(dst);
	for(int i=result.size()-1; i>=0; i--) cout << result[i] << " ";
}

void dijkstraAlgorithm(int src, int dst){
	priority_queue<edge*, vector<edge*>, Compare> pq;
	d[src] = 0;
	pq.push(new edge(src, src, 0));
	while(pq.size() > 0){
    edge* e1 = pq.top(); pq.pop();
    int u = e1->u;
		visited[u] = true;
		if(visited[dst] == true) break;
		string blockName = (nodeid_cordinates[u])->block_x_y;
		readFromFileHelperFunction(blockName);
		for(int i=0; i<adj[u].size(); i++){
			pair<int, long double> temp = adj[u][i];
			if(d[u]+temp.second < d[temp.first]){
				d[temp.first] = d[u]+temp.second;
				parent[temp.first] = u;
				pq.push(new edge(temp.first, temp.first, d[temp.first]));
			}
		}
  }
	if(d[dst] == INT_MAX){
		cout << "Distance: INF";
	} else{
		cout << setprecision (numeric_limits<double>::digits10 + 1) << "Distance: " << d[dst] << " ";
		cout << " \tPath: ";
		printPath(dst);
		cout << "\tBlocks read: " << blockCount;
	}
  cout << "\n";
}

void shortestPathHelperFunction(){
	int src, dst;
	cout << "Enter Source and Destination node:  ";
	cin >> src >> dst;
	if(nodeid_cordinates.find(src) == nodeid_cordinates.end() || nodeid_cordinates.find(dst) == nodeid_cordinates.end()){
		cout << "Source or Destination not present\n";
		return;
	}
	string blockName = (nodeid_cordinates[src])->block_x_y;
	readFromFileHelperFunction(blockName);
	dijkstraAlgorithm(src, dst);
	for(map<int, bool>::iterator it = visited.begin(); it != visited.end(); it++) visited[it->first] = false;
	for(map<int, long double>::iterator it = d.begin(); it != d.end(); it++) d[it->first] = INT_MAX;
	for(map<int, int>::iterator it = parent.begin(); it != parent.end(); it++) parent[it->first] = -1;
	adj.clear();
	filesRead.clear();
	blockCount = 0;
}

int main(){
	int option;
	cout << "Single shortest path algorithm implementation - Dijkstra Algorithm\n";
	cout << "1. Create Partition\t 2. Find the shortest path\t 3. Exit\n";
	cout << "Choose an option:  ";
	cin >> option;
	while(true){
		switch (option) {
			case /* value */1: partitionHelperFunction(); break;
			case /* value */2: shortestPathHelperFunction(); break;
			default: return 0;
		}
		cout << "Choose an option:  ";
		cin >> option;
	}
	return 0;
}
