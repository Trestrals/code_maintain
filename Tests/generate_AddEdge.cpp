#include <bits/stdc++.h>
using namespace std;



void writeGraph(vector<pair<int, int> > edges, int n, string fname) {
    FILE* f = freopen(fname.c_str(), "w", stdout);
    if(n != 0) printf("%d %d\n", n, edges.size());
    for(auto e : edges) {
        printf("%d %d\n", e.first, e.second);
    }
    fclose(f);
}

int numN[10000000];
unordered_map<int, bool> choose;

int main(int argc, char *argv[]){

    if(argc != 3) exit(-1);

    srand(time(0));
    string name = argv[1];

    string path = "../datasets/_" + name + ".txt";

    freopen(path.c_str(), "r", stdin);

    int n = 0;
    int m = 0;

    scanf("%d %d", &n, &m);

    vector<pair<int, int> > edge;

    for(int i = 1; i <= m; i++){
        int a, b;
        scanf("%d %d", &a, &b);
        edge.push_back(make_pair(a, b));
        numN[a]++;
        numN[b]++;
    }

    vector<pair<int, int> > edgePre;
    vector<pair<int, int> > edgeAdd;

    int chs = stoi(argv[2]);

    for(int i = 1; i <= chs; i++){
        int x = rand() % m;
        while(numN[edge[x].first] == 1 || numN[edge[x].second] == 1 || choose[x]) {
            x = rand() % m;
        }
        numN[edge[x].first]--;
        numN[edge[x].second]--;
        choose[x] = true;
    }

    for(int i = 0; i < m; i++){
        if(choose[i]) edgeAdd.push_back(edge[i]);
        else edgePre.push_back(edge[i]);
    }

    random_shuffle(edgeAdd.begin(), edgeAdd.end());

    string path1 = "../datasets/_" + name + "2.txt";
    writeGraph(edgePre, n, path1);
    string path2 = "../datasets/_" + name + "2_addGraph.txt";
    writeGraph(edgeAdd, 0, path2);


    return 0;

}