#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <sstream>
#include <limits>

using namespace std;

/*

Forma optima de cablear con fibra optica,
de forma que se pueda compartir informacion entre dos 
colonias cuales quiera.

Por medio de Kruskal -->

*/

// Estructura para representar una arista
struct Edge {
    int u, v;
    int weight;
    bool operator<(const Edge& other) const {
        return weight < other.weight;
    }
};

// Función para encontrar el representante de un nodo (para el Union-Find)
int findSet(int node, vector<int>& parent) {
    if (parent[node] != node)
        parent[node] = findSet(parent[node], parent); // Path compression
    return parent[node];
}

// Función para unir dos conjuntos (para el Union-Find)
void unionSets(int u, int v, vector<int>& parent, vector<int>& rank) {
    int rootU = findSet(u, parent);
    int rootV = findSet(v, parent);
    if (rootU != rootV) {
        if (rank[rootU] > rank[rootV]) {
            parent[rootV] = rootU;
        } else if (rank[rootU] < rank[rootV]) {
            parent[rootU] = rootV;
        } else {
            parent[rootV] = rootU;
            rank[rootU]++;
        }
    }
}

// Función para calcular el MST usando el algoritmo de Kruskal
vector<Edge> kruskalMST(int n, vector<Edge>& edges) {
    sort(edges.begin(), edges.end()); // Ordenar las aristas por peso
    vector<int> parent(n), rank(n, 0);
    vector<Edge> mst;

    // Inicializar el Union-Find
    for (int i = 0; i < n; ++i)
        parent[i] = i;

    // Procesar las aristas en orden de peso
    for (const auto& edge : edges) {
        if (findSet(edge.u, parent) != findSet(edge.v, parent)) {
            mst.push_back(edge);
            unionSets(edge.u, edge.v, parent, rank);
        }
    }
    return mst;
}

// Función para leer la matriz desde un archivo
vector<Edge> readGraph(const string& filename, int& n) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error al abrir el archivo de entrada." << endl;
        exit(1);
    }

    file >> n; // Leer el número de nodos
    vector<Edge> edges;

    // Leer la matriz y convertirla en una lista de aristas
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            int weight;
            file >> weight;
            if (i < j && weight > 0) { // Solo tomar la parte superior de la matriz
                edges.push_back({i, j, weight});
            }
        }
    }
    return edges;
}


/*

Ruta mas corta posible, que visita cada colonia
exactamente una vez y al finalizar regresa a 
la colonia origen.

Por medio de branch and bound -->

*/


int main(int argc, char* argv[]) {
    if (argc != 2) {
        cerr << "Uso: " << argv[0] << " <archivo_entrada>" << endl;
        return 1;
    }

    string filename = argv[1];
    int n; // Número de nodos
    vector<Edge> edges = readGraph(filename, n);

    // Calcular el MST
    vector<Edge> mst = kruskalMST(n, edges);

    // Imprimir las aristas del MST
    cout << "Forma de cablear las colonias con fibra (MST):" << endl;
    for (const auto& edge : mst) {
        cout << "(" << edge.u << ", " << edge.v << ") -> " << edge.weight << " km" << endl;
    }

    return 0;
}
