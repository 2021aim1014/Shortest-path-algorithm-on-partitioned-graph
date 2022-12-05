#include<bits/stdc++.h>
#include <fstream>
using namespace std;

// #define NODES "testNodes.txt"
// #define EDGES "testEdges.txt"

#define NODES "sample_nodes.txt"
#define EDGES "sample_edges.txt"

// #define NODES "nodes.txt"
// #define EDGES "edges.txt"

struct edge{
	int u, v;
	long double length;
	edge(int u, int v, long double length){
		this->u = u;
		this->v = v;
		this->length = length;
	}
};

class Compare{
public:
  int operator() (edge* v1, edge* v2){
    return v1->length > v2->length;
  }
};

void printPath(int dst, map<int, int> &parent){
	vector<int> result;
	while(parent[dst] != -1){
		result.push_back(dst);
		dst = parent[dst];
	}
	result.push_back(dst);
	for(int i=result.size()-1; i>=0; i--) cout << result[i] << " ";
}

int dijkstra(int src, int dst, map<int, long double> &d, map<int, vector<pair<int, long double> > > &adj, map<int, bool> &visited, map<int, int> &parent){
	d[src] = 0;
	priority_queue<edge*, vector<edge*>, Compare> pq;
	pq.push(new edge(src, src, 0));
	while(pq.size() > 0){
		edge* e1 = pq.top(); pq.pop();
		int u = e1->u;
		visited[u] = true;
		if(visited[dst] == true) break;
		for(int i=0; i<adj[u].size(); i++){
			pair<int, long double> temp = adj[u][i];
			if(d[u]+temp.second < d[temp.first]){
				d[temp.first] = d[u]+temp.second;
				parent[temp.first] = u;
				pq.push(new edge(temp.first, temp.first, d[temp.first]));
			}
		}
	}
	if(d[dst] == DBL_MAX){
		cout << "Distance: INF";
	} else{
		cout << setprecision (numeric_limits<double>::digits10 + 1) << "Distance: " << d[dst] << " ";
		cout << "  Path: ";
		printPath(dst, parent);
	}
  cout << "\n";
	return 0;
}

int main(){

  map<int, long double> d;
  map<int, vector<pair<int, long double> > > adj;
	map<int, bool> visited;
	map<int, int> parent;
	ifstream ifp1(NODES);
	if(!ifp1) return -1;
	while(!ifp1.eof()){
		int node;
		long double x, y;
		ifp1 >> node >> x >> y;
		//c
		visited[node] = false;
		d[node] = DBL_MAX;
		parent[node] = -1;
	}
	ifp1.close();

	ifstream ifp2(EDGES);
	if(!ifp2) return -1;
	set<string> edgesDuplicates;
	while(!ifp2.eof()){
		int u, v;
		long double w;
		ifp2 >> u >> v >> w;
		string temp1 = to_string(u)+"_"+to_string(v);
		string temp2 = to_string(v)+"_"+to_string(u);
		if(edgesDuplicates.find(temp1) != edgesDuplicates.end()) continue;
		edgesDuplicates.insert(temp1);
		edgesDuplicates.insert(temp2);
    pair<int, long double> vw, uw;
    vw.first = v; vw.second = w; uw.first = u; uw.second = w;
    adj[u].push_back(vw);
    adj[v].push_back(uw);
	}
	ifp2.close();
	while(true){
		int src, dst;
		cin >> src >> dst;
	  dijkstra(src, dst, d, adj, visited, parent);
		for(map<int, bool>::iterator it = visited.begin(); it != visited.end(); it++) visited[it->first] = false;
		for(map<int, long double>::iterator it = d.begin(); it != d.end(); it++) d[it->first] = DBL_MAX;
		for(map<int, int>::iterator it = parent.begin(); it != parent.end(); it++) parent[it->first] = -1;
	}
  return 0;
}
