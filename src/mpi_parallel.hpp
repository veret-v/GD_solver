#pragma once

#include <mpi.h>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>
#include <limits>
#include <utility>

#include "mesh.hpp"

namespace cfd {

struct HaloComm {
    int remote_rank;
    std::vector<int> send_ids; // Индексы наших ячеек для отправки
    std::vector<int> recv_ids; // Индексы ghost-ячеек для приема
};

struct Halo {
    std::vector<HaloComm> comms;
    MPI_Datatype cell_mpi_type;
};

// Плоская структура только с POD-данными для MPI
struct CellDataMPI {
    double U[4];
};

inline MPI_Datatype create_cell_mpi_type() {
    CellDataMPI dummy{};
    MPI_Aint base, dispU;
    MPI_Get_address(&dummy, &base);
    MPI_Get_address(&dummy.U, &dispU);
    
    MPI_Aint disp[1] = {dispU - base};
    int blocks[1] = {4};
    MPI_Datatype types[1] = {MPI_DOUBLE};

    MPI_Datatype raw_type, final_type;
    MPI_Type_create_struct(1, blocks, disp, types, &raw_type);
    
    // Обязательный resized для защиты от padding'а
    MPI_Type_create_resized(raw_type, 0, sizeof(CellDataMPI), &final_type);
    MPI_Type_commit(&final_type);
    MPI_Type_free(&raw_type);
    return final_type;
}

inline bool less_by_axis(const std::vector<Cell>& cells, bool cut_x, int a, int b) {
    const double ca = cut_x ? cells[a].cx : cells[a].cy;
    const double cb = cut_x ? cells[b].cx : cells[b].cy;
    if (ca != cb) return ca < cb;

    const double oa = cut_x ? cells[a].cy : cells[a].cx;
    const double ob = cut_x ? cells[b].cy : cells[b].cx;
    if (oa != ob) return oa < ob;

    return a < b;
}

inline int rcb_split_index(int count, int K_left, int K) {
    int split = static_cast<int>(std::round(static_cast<double>(count) * K_left / K));
    if (count > 1) split = std::clamp(split, 1, count - 1);
    return split;
}

inline void partition_rcb_range(std::vector<int>& ids, int first, int last, int K,
                                int off, const std::vector<Cell>& cells,
                                std::vector<int>& part) {
    const int count = last - first;
    if (K <= 1 || count <= 1) {
        for (int i = first; i < last; ++i) part[ids[i]] = off;
        return;
    }

    const int active_K = std::min(K, count);
    const int K_left = active_K / 2;
    const int K_right = active_K - K_left;

    double xmin = std::numeric_limits<double>::max();
    double xmax = std::numeric_limits<double>::lowest();
    double ymin = std::numeric_limits<double>::max();
    double ymax = std::numeric_limits<double>::lowest();

    for (int i = first; i < last; ++i) {
        const Cell& c = cells[ids[i]];
        xmin = std::min(xmin, c.cx);
        xmax = std::max(xmax, c.cx);
        ymin = std::min(ymin, c.cy);
        ymax = std::max(ymax, c.cy);
    }
    
    // Выбираем самую длинную ось
    bool cut_x = (xmax - xmin) >= (ymax - ymin);
    
    // Сортируем по центроидам
    std::sort(ids.begin() + first, ids.begin() + last, [&](int a, int b) {
        return less_by_axis(cells, cut_x, a, b);
    });
    
    // Режем пропорционально числу процессов
    const int split = first + rcb_split_index(count, K_left, active_K);

    partition_rcb_range(ids, first, split, K_left, off, cells, part);
    partition_rcb_range(ids, split, last, K_right, off + K_left, cells, part);
}

// CHECK: DECOMPOSITION
// 1. Алгоритм RCB (Рекурсивная координатная бисекция)
inline void partition_rcb(std::vector<int>& ids, int K, int off, const std::vector<Cell>& cells, std::vector<int>& part) {
    partition_rcb_range(ids, 0, static_cast<int>(ids.size()), K, off, cells, part);
}

// 2. Сборка локальной сетки и определение MPI-границ
inline Mesh build_local_mesh(const Mesh& global, const std::vector<int>& part, int rank, Halo& halo) {
    Mesh m;
    m.nodes = global.nodes; // Для простоты копируем все узлы

    std::vector<int> g2l(global.nc(), -1);
    int n_owned = 0;

    // Сначала добавляем только "свои" ячейки
    for(int i = 0; i < global.nc(); ++i) {
        if(part[i] == rank) {
            g2l[i] = n_owned++;
            m.cells.push_back(global.cells[i]);
        }
    }
    m.n_owned = n_owned;

    int n_ghosts = 0;
    struct PendingHaloComm {
        std::vector<std::pair<int, int>> send_ids; // global id, local id
        std::vector<std::pair<int, int>> recv_ids; // global id, local id
    };

    std::map<int, PendingHaloComm> cmap;

    // Проходим по глобальным граням и собираем локальные + ghost
    for(const Face& gf : global.faces) {
        int pL = part[gf.left];
        int pR = gf.right >= 0 ? part[gf.right] : -1;

        if (pL != rank && pR != rank) continue; // Грань вообще не наша

        Face lf = gf;
        if (pL == rank && pR == rank) {
            // Внутренняя грань
            lf.left = g2l[gf.left];
            lf.right = g2l[gf.right];
            lf.remote_rank = -1;
            lf.bc = Face::BC::Interior;
        } else if (pL == rank && pR != rank) {
            lf.left = g2l[gf.left];
            if (pR >= 0) { 
                // MPI граница: мы слева, сосед справа
                if (g2l[gf.right] == -1) {
                    g2l[gf.right] = n_owned + n_ghosts++;
                    m.cells.push_back(global.cells[gf.right]); // Добавляем ghost
                }
                lf.right = g2l[gf.right];
                lf.remote_rank = pR;
                lf.bc = Face::BC::MPIBound;
                cmap[pR].send_ids.push_back({gf.left, lf.left});
                cmap[pR].recv_ids.push_back({gf.right, lf.right});
            }
        } else if (pR == rank && pL != rank) {
            // MPI граница: мы справа, сосед слева. Переворачиваем нормаль, чтобы своя ячейка всегда была left!
            lf.left = g2l[gf.right];
            if (g2l[gf.left] == -1) {
                g2l[gf.left] = n_owned + n_ghosts++;
                m.cells.push_back(global.cells[gf.left]); // Добавляем ghost
            }
            lf.right = g2l[gf.left];
            lf.nx = -lf.nx;
            lf.ny = -lf.ny;
            lf.remote_rank = pL;
            lf.bc = Face::BC::MPIBound;
            cmap[pL].send_ids.push_back({gf.right, lf.left});
            cmap[pL].recv_ids.push_back({gf.left, lf.right});
        }
        
        m.faces.push_back(lf);
    }

    // Перестраиваем face_ids для локальных ячеек
    for(auto& c : m.cells) c.face_ids.clear();
    for(int fid = 0; fid < static_cast<int>(m.faces.size()); ++fid) {
        m.cells[m.faces[fid].left].face_ids.push_back(fid);
        if (m.faces[fid].right >= 0) {
            m.cells[m.faces[fid].right].face_ids.push_back(fid);
        }
    }

    auto ordered_local_ids = [](std::vector<std::pair<int, int>>& ids) {
        std::sort(ids.begin(), ids.end(), [](const auto& a, const auto& b) {
            return a.first < b.first;
        });
        ids.erase(std::unique(ids.begin(), ids.end(), [](const auto& a, const auto& b) {
            return a.first == b.first;
        }), ids.end());

        std::vector<int> local_ids;
        local_ids.reserve(ids.size());
        for (const auto& id : ids) local_ids.push_back(id.second);
        return local_ids;
    };

    // Завершаем сборку Halo
    for(auto& [r, hc] : cmap) {
        HaloComm comm;
        comm.remote_rank = r;
        // Удаляем дубликаты ячеек, т.к. процессы могут граничить по нескольким граням
        comm.send_ids = ordered_local_ids(hc.send_ids);
        comm.recv_ids = ordered_local_ids(hc.recv_ids);
        halo.comms.push_back(comm);
    }
    
    halo.cell_mpi_type = create_cell_mpi_type();
    return m;
}

// CHECK: HALO_EXCHANGE
// 3. Неблокирующий обмен Halo
inline void exchange_halo(Mesh& m, Halo& halo, MPI_Comm comm) {
    std::vector<MPI_Request> reqs;
    std::vector<std::vector<CellDataMPI>> send_bufs(halo.comms.size());
    std::vector<std::vector<CellDataMPI>> recv_bufs(halo.comms.size());

    // 1. Постим все Irecv
    int idx = 0;
    for(auto& hc : halo.comms) {
        recv_bufs[idx].resize(hc.recv_ids.size());
        MPI_Request req;
        MPI_Irecv(recv_bufs[idx].data(), recv_bufs[idx].size(), halo.cell_mpi_type,
                  hc.remote_rank, 0, comm, &req);
        reqs.push_back(req);
        idx++;
    }

    // 2. Пакуем буферы и постим Isend
    idx = 0;
    for(auto& hc : halo.comms) {
        send_bufs[idx].resize(hc.send_ids.size());
        for(size_t i = 0; i < hc.send_ids.size(); ++i) {
            for(int k = 0; k < 4; ++k) send_bufs[idx][i].U[k] = m.cells[hc.send_ids[i]].U[k];
        }
        MPI_Request req;
        MPI_Isend(send_bufs[idx].data(), send_bufs[idx].size(), halo.cell_mpi_type,
                  hc.remote_rank, 0, comm, &req);
        reqs.push_back(req);
        idx++;
    }

    // 3. Ждем завершения всех операций
    MPI_Waitall(reqs.size(), reqs.data(), MPI_STATUSES_IGNORE);

    // 4. Распаковываем данные в Ghost-ячейки
    idx = 0;
    for(auto& hc : halo.comms) {
        for(size_t i = 0; i < hc.recv_ids.size(); ++i) {
            for(int k = 0; k < 4; ++k) m.cells[hc.recv_ids[i]].U[k] = recv_bufs[idx][i].U[k];
        }
        idx++;
    }
}

} // namespace cfd
