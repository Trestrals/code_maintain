#include <bits/stdc++.h>
using namespace std;

void writeGraph(vector<pair<int, int> > edges, int n, string fname) {
    // puts("IN");
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
    string pstr = argv[2];
    double p = stod(pstr) / 100;

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

    vector<int> edgechs;
    vector<int> vertex;
    edgechs.resize(m, -1);
    vertex.resize(n, -1);

    int chs = m * p;

    for(int i = 0; i < chs; i++){
        int x = rand() % m;
        while(edgechs[x] != -1) {
            x = rand() % m;
        }
        edgechs[x] = 1;
        vertex[edge[x].first] = 1;
        vertex[edge[x].second] = 1;
    }

    int cntn = 0;
    for (int i = 0; i < n; i++) {
        if (vertex[i] == 1) {
            vertex[i] = cntn;
            cntn++;
        }
    }

    for(int i = 0; i < m; i++){

        int a = edge[i].first;
        int b = edge[i].second;
        if(edgechs[i] == 1) {
            edgePre.push_back({vertex[a], vertex[b]});
        }
    }

    // puts("HAA");
    cout << cntn << " " << edgePre.size() << endl;

    // random_shuffle(edgeAdd.begin(), edgeAdd.end());
    // string path1 = name + "aa.txt";

   string path1 = "../datasets/_" + name + "_" + pstr + "E.txt";
    writeGraph(edgePre, cntn, path1);


    return 0;

}