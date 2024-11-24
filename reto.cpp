#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include <queue>
#include <fstream>
#include <climits>
#include <algorithm>
#include <limits>

using namespace std;

/*
Forma optima de cablear con fibra optica,
de forma que se pueda compartir informacion entre dos 
colonias cuales quiera.

Por medio de PRIM -->
*/
class PrimMST {
private:
    int V; // Numero de vertices
    vector<vector<int>> graph; // Matriz de adyacencia

public:
    PrimMST(int vertices) : V(vertices) {
        graph.resize(V, vector<int>(V));
    }

    void setGraph(const vector<vector<int>>& g) { graph = g; }

    // Funcion para encontrar el vertice con la distancia minima
    int minKey(vector<int>& key, vector<bool>& mstSet) {
        int min = INT_MAX, min_index;
        for (int v = 0; v < V; v++) {
            if (!mstSet[v] && key[v] < min) {
                min = key[v];
                min_index = v;
            }
        }
        return min_index;
    }

    // Imprimir el MST
    void printMST(vector<int>& parent) {
        cout << "\n1. Forma de cablear las colonias:\n" << endl;
        int totalCost = 0;
        for (int i = 1; i < V; i++) {
            char origen = 'A' + parent[i];
            char destino = 'A' + i;
            cout << "(" << origen << "," << destino << ")" << endl;
            totalCost += graph[i][parent[i]];
        }
        cout << "\nCosto total del cableado: " << totalCost << " km" << endl;
    }

    // PRIM
    void findMST() {
        vector<int> parent(V);     // Array para almacenar el MST
        vector<int> key(V);        // Valores clave para elegir el peso minimo
        vector<bool> mstSet(V);    // Para representar el conjunto de vertices incluidos

        // Inicializar como infinito
        for (int i = 0; i < V; i++) {
            key[i] = INT_MAX;
            mstSet[i] = false;
        }

        key[0] = 0; // Primer nodo como raiz
        parent[0] = -1; 

        for (int count = 0; count < V - 1; count++) {
            int u = minKey(key, mstSet);
            mstSet[u] = true;

            for (int v = 0; v < V; v++) {
                if (graph[u][v] && !mstSet[v] && graph[u][v] < key[v]) {
                    parent[v] = u;
                    key[v] = graph[u][v];
                }
            }
        }
        printMST(parent);
    }
};


/*
Ruta mas corta posible, que visita cada colonia
exactamente una vez y al finalizar regresa a 
la colonia origen.

Por medio de HELD-KARP -->
*/
class HeldKarpTSP {
private:
    int n;
    vector<vector<int>> dist;
    map<pair<int, int>, int> memo;

    int dp(int pos, int mask) {
        if (mask == ((1 << n) - 1)) {
            return dist[pos][0]; // Retorno a la ciudad inicial
        }

        pair<int, int> state = {pos, mask};
        if (memo.find(state) != memo.end()) {
            return memo[state];
        }

        int ans = INT_MAX;
        for (int city = 0; city < n; city++) {
            if (!(mask & (1 << city))) {
                int newAns = dist[pos][city] + dp(city, mask | (1 << city));
                ans = min(ans, newAns);
            }
        }

        return memo[state] = ans;
    }

    void findPath(int pos, int mask, vector<int>& path) {
        if (mask == ((1 << n) - 1)) {
            path.push_back(0);
            return;
        }

        int minCost = INT_MAX, nextCity = -1;

        for (int city = 0; city < n; city++) {
            if (!(mask & (1 << city))) {
                int cost = dist[pos][city] + dp(city, mask | (1 << city));
                if (cost < minCost) {
                    minCost = cost;
                    nextCity = city;
                }
            }
        }

        path.push_back(nextCity);
        findPath(nextCity, mask | (1 << nextCity), path);
    }

public:
    void setDistances(const vector<vector<int>>& distances) {
        dist = distances;
        n = dist.size();
    }

    void solve() {
        memo.clear();
        vector<int> path = {0}; // Comenzamos en ciudad 0 - A
        int minCost = dp(0, 1); // Ciudad 0 con mascara inicial 1
        findPath(0, 1, path); // Encontrar camino

        cout << "\n----------------------------------------" << endl; 
        cout << "2. Ruta a seguir por el personal:\n" << endl;
        for (size_t i = 0; i < path.size(); i++) {
            cout << char('A' + path[i]);
            if (i < path.size() - 1) cout << " -> ";
        }
        cout << endl;
        cout << "\nCosto total del recorrido: " << minCost << " km" << endl;
    }
};


/*
Conocer el flujo maximo de informacion entre 
una colonia i y una colonia j, o nodo inicial 
al nodo final.

Por medio de EDMONDS-KARP -->
*/
class EdmondsKarp {
private:
    int V; //Numero de vertices
    vector<vector<int>> capacity; //Matriz de capacidades de flujo

    bool bfs(vector<vector<int>>& residual, vector<int>& parent, int s, int t) {
        vector<bool> visited(V, false);
        queue<int> q;
        q.push(s);
        visited[s] = true;
        parent[s] = -1;

        while (!q.empty()) {
            int u = q.front();
            q.pop();

            for (int v = 0; v < V; v++) {
                if (!visited[v] && residual[u][v] > 0) {
                    q.push(v);
                    parent[v] = u;
                    visited[v] = true;
                    if (v == t) return true;
                }
            }
        }
        return false;
    }

public:
    void setCapacities(const vector<vector<int>>& capacities) {
        capacity = capacities;
        V = capacity.size();
    }

    int maxFlow(int source, int sink) {
        vector<vector<int>> residual = capacity;
        vector<int> parent(V);
        int max_flow = 0;

        while (bfs(residual, parent, source, sink)) {
            int path_flow = INT_MAX;
            for (int v = sink; v != source; v = parent[v]) {
                int u = parent[v];
                path_flow = min(path_flow, residual[u][v]);
            }
            for (int v = sink; v != source; v = parent[v]) {
                int u = parent[v];
                residual[u][v] -= path_flow;
                residual[v][u] += path_flow;
            }
            max_flow += path_flow;
        }
        return max_flow;
    }
};


/*
Dada una nueva contratacion del servicio, 
cual es la central mas cercana geograficamente 
a esa nueva contratacion.

Por medio de VORONOI -->
*/
struct Point {
    double x, y;
    Point(double _x = 0, double _y = 0) : x(_x), y(_y) {}
};

class VoronoiDiagram {
private:
    vector<Point> centrals;

    double calculateDistance(const Point& p1, const Point& p2) {
        double dx = p1.x - p2.x;
        double dy = p1.y - p2.y;
        return sqrt(dx * dx + dy * dy);
    }

public:
    void setCentrals(const vector<Point>& points) { centrals = points; }

    void findNearestCentral(const Point& newLocation) {
        int nearestIdx = 0;
        double minDist = calculateDistance(newLocation, centrals[0]);

        for (size_t i = 1; i < centrals.size(); i++) {
            double dist = calculateDistance(newLocation, centrals[i]);
            if (dist < minDist) {
                minDist = dist;
                nearestIdx = i;
            }
        }

        cout << "\n----------------------------------------" << endl;
        cout << "4. Central mas cercana:\n" << endl;
        cout << "Coordenadas ingresadas: (" << newLocation.x << ", " << newLocation.y << ")" << endl;
        cout << "Central mas cercana: " << (char)('A' + nearestIdx) << " en (" << centrals[nearestIdx].x << ", " << centrals[nearestIdx].y << ")" << endl;
        cout << "Distancia: " << minDist << "\n" << endl;
    }
};


// LEER EL ARCHIVO //
void readInputFile(const string& filename, vector<vector<int>>& graph, vector<vector<int>>& capacities, vector<Point>& centrals) {
    ifstream file(filename);
    if (!file.is_open()) {
        cout << "Error al abrir el archivo" << endl;
        exit(1);
    }

    int n;
    file >> n;

    graph.resize(n, vector<int>(n));
    capacities.resize(n, vector<int>(n));

    // Leer matriz de adyacencia Prim y Held-Karp
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            file >> graph[i][j];
        }
    }

    // Leer matriz de capacidades Edmonds-Karp
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            file >> capacities[i][j];
        }
    }

    // Leer puntos para Voronoi
    char dummy;
    double x, y;
    while (file >> dummy >> x >> dummy >> y >> dummy) {
        centrals.push_back(Point(x, y));
    }

    file.close();
}

// --- MAIN ---
int main(int argc, char* argv[]) {
    if (argc != 6) {
        cout << "Ejemplo de uso: " << " <archivo.txt> nodoInicial(0 a n) nodoFinal(0 a n) coordenadaX coordenadaY" << endl;
        return 1;
    }

    string filename = argv[1];
    int nodoinicial = atoi(argv[2]);
    int nodofinal = atoi(argv[3]);
    double coordX = atof(argv[4]);
    double coordY = atof(argv[5]);

    vector<vector<int>> graph, capacities;
    vector<Point> centrals;

    readInputFile(filename, graph, capacities, centrals);

    // Problema 1
    PrimMST prim(graph.size());
    prim.setGraph(graph);
    prim.findMST();

    // Problema 2
    HeldKarpTSP tsp;
    tsp.setDistances(graph);
    tsp.solve();

    // Problema 3
    EdmondsKarp ek;
    ek.setCapacities(capacities);
    cout << "\n----------------------------------------" << endl;
    cout << "3. Flujo maximo de informacion:\n" << endl;
    cout << "Desde el nodo " << nodoinicial << " hasta el nodo " << nodofinal << ": ";
    cout << ek.maxFlow(nodoinicial, nodofinal) << endl;

    // Problema 4
    VoronoiDiagram vd;
    vd.setCentrals(centrals);
    vd.findNearestCentral(Point(coordX, coordY));

    return 0;
}
