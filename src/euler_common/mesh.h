#pragma once
#include <vector>
#include <array>

using Vec4 = std::array<double, 4>;

struct Node {
    double x, y;
};

struct Face {
    int left, right; // индексы ячеек слева и справа
    double nx, ny;   // вектор нормали
    double length;   // длина грани
    double mx, my;   // центр грани
    int bc_type;     // код граничного условия (аналог ваших bc_code)
    
    bool is_boundary() const { return right < 0; }
};

struct Cell {
    std::vector<int> face_ids;
    std::vector<int> node_ids;
    double vol; 
    double cx, cy; 
    Vec4 U;    
    Vec4 res; 
    int partition; // Для METIS
};

struct CommTarget {
    int rank;
    std::vector<int> send_indices; 
    std::vector<int> recv_indices; 
};

struct Mesh {
    std::vector<Node> nodes;
    std::vector<Cell> cells;
    std::vector<Face> faces;
    std::vector<CommTarget> comm_targets; 
    
    int nc() const { return cells.size(); }
    int nf() const { return faces.size(); }
    void zero_res() {
        for (auto& c : cells) c.res.fill(0.0);
    }
};