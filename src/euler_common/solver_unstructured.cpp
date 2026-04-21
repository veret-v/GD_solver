#include "mesh.h"
#include <mpi.h>
#include <cmath>
#include <algorithm>
#include <iostream>

// Используем вашу функцию HLLC с нормалями, которую я адаптировал ранее
extern std::vector<double> hllc_flux_unstructured(const Vec4& UL, const Vec4& UR, 
                                                  double nx, double ny, double g);

// Функция генерации "призрачного" состояния для границ
Vec4 ghost_state(const Vec4& UL, double nx, double ny, int bc_type, double g) {
    // По умолчанию жесткая стенка (Slip wall)
    double rho = UL[0];
    double u = UL[1] / rho;
    double v = UL[2] / rho;
    
    double un = u * nx + v * ny;
    
    Vec4 ghost = UL;
    ghost[1] = rho * (u - 2.0 * un * nx);
    ghost[2] = rho * (v - 2.0 * un * ny);
    
    return ghost;
}

void compute_residuals(Mesh& m, double g, int rank) {
    m.zero_res();
    
    for (const auto& f : m.faces) {
        // Пропускаем грани, которые не относятся к нашему MPI-процессу
        if (m.cells[f.left].partition != rank && 
           (f.right >= 0 && m.cells[f.right].partition != rank)) {
            continue; 
        }

        Vec4 UL = m.cells[f.left].U;
        Vec4 UR;

        if (f.is_boundary()) {
            UR = ghost_state(UL, f.nx, f.ny, f.bc_type, g); 
        } else {
            UR = m.cells[f.right].U;
        }

        // Расчет потока через грань
        std::vector<double> flux = hllc_flux_unstructured(UL, UR, f.nx, f.ny, g);

        // Распределяем потоки по ячейкам
        for(int k = 0; k < 4; ++k) {
            if (m.cells[f.left].partition == rank) {
                m.cells[f.left].res[k] -= flux[k] * f.length;
            }
            if (f.right >= 0 && m.cells[f.right].partition == rank) {
                m.cells[f.right].res[k] += flux[k] * f.length;
            }
        }
    }
}

// MPI Обмен гало-ячейками (теневыми)
void exchange_unstructured_halos(Mesh& m) {
    std::vector<MPI_Request> requests;
    std::vector<std::vector<double>> send_bufs(m.comm_targets.size());
    std::vector<std::vector<double>> recv_bufs(m.comm_targets.size());

    // Инициируем асинхронные обмены
    for (size_t i = 0; i < m.comm_targets.size(); ++i) {
        const auto& target = m.comm_targets[i];
        send_bufs[i].resize(target.send_indices.size() * 4);
        recv_bufs[i].resize(target.recv_indices.size() * 4);

        for (size_t j = 0; j < target.send_indices.size(); ++j) {
            int cell_idx = target.send_indices[j];
            for(int k=0; k<4; ++k) send_bufs[i][j*4 + k] = m.cells[cell_idx].U[k];
        }

        MPI_Request req_recv, req_send;
        MPI_Irecv(recv_bufs[i].data(), recv_bufs[i].size(), MPI_DOUBLE, 
                  target.rank, 0, MPI_COMM_WORLD, &req_recv);
        MPI_Isend(send_bufs[i].data(), send_bufs[i].size(), MPI_DOUBLE, 
                  target.rank, 0, MPI_COMM_WORLD, &req_send);
                  
        requests.push_back(req_recv);
        requests.push_back(req_send);
    }

    MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);

    // Распаковываем полученные данные
    for (size_t i = 0; i < m.comm_targets.size(); ++i) {
        const auto& target = m.comm_targets[i];
        for (size_t j = 0; j < target.recv_indices.size(); ++j) {
            int cell_idx = target.recv_indices[j];
            for(int k=0; k<4; ++k) m.cells[cell_idx].U[k] = recv_bufs[i][j*4 + k];
        }
    }
}

// Главный цикл для неструктурированной сетки
void solve_unstructured(int rank, int size, double tmax, double cfl, double g, 
                        const std::string& mesh_file) {
    Mesh mesh;
    
    // В реальном коде здесь вызывается функция чтения .msh файла через Gmsh API
    // load_gmsh_mesh(mesh_file, mesh);
    
    // Если запускаем в параллели, вызываем METIS
    if (size > 1 && rank == 0) {
        // partition_mesh_metis(mesh, size);
        // distribute_mesh_via_mpi(mesh);
    }

    double curr_time = 0.0;
    int step = 0;

    while (curr_time < tmax) {
        // 1. Обмен границами
        if (size > 1) {
            exchange_unstructured_halos(mesh);
        }

        // 2. Расчет шага по времени (упрощенно)
        double local_dt = 1e10; // Вычисляется как V_i / (L_f * (|u_n| + a))
        double dt;
        MPI_Allreduce(&local_dt, &dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        if (curr_time + dt > tmax) dt = tmax - curr_time;

        // 3. Вычисление потоков и невязок
        compute_residuals(mesh, g, rank);

        // 4. Явный шаг по времени
        for (auto& c : mesh.cells) {
            if (c.partition == rank) {
                for (int k = 0; k < 4; ++k) {
                    c.U[k] += (dt / c.vol) * c.res[k];
                }
            }
        }

        curr_time += dt;
        step++;
        
        if (rank == 0 && step % 100 == 0) {
            std::cout << "Unstructured Step: " << step << " Time: " << curr_time << std::endl;
        }
    }
}