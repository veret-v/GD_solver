#include "script.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <regex>

void parseAndInitialize(std::vector<std::vector<std::vector<double>>>& u, 
                        const std::string& filename, 
                        int offset_x, 
                        int offset_y, 
                        double gamma) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "ERROR: Cannot open " << filename << std::endl;
        return;
    }

    std::stringstream buffer;
    buffer << file.rdbuf();
    std::string content = buffer.str();

    std::regex block_pattern(R"(BLOCK_\d+([\s\S]*?)(?=BLOCK_\d+|$))");
    // Обновленная регулярка: поддерживает и "=", и ":"
    std::regex kv_pattern(R"(([\w_]+)\s*[:=]\s*([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?))");

    auto blocks_begin = std::sregex_iterator(content.begin(), content.end(), block_pattern);
    auto blocks_end = std::sregex_iterator();

    for (auto i = blocks_begin; i != blocks_end; ++i) {
        std::string block_body = (*i)[1].str();
        int xb = 0, xe = 0, yb = 0, ye = 0;
        double rho = 0, u_vel = 0, v_vel = 0, p = 0;

        auto kv_begin = std::sregex_iterator(block_body.begin(), block_body.end(), kv_pattern);
        for (auto j = kv_begin; j != std::sregex_iterator(); ++j) {
            std::string key = (*j)[1].str();
            double val = std::stod((*j)[2].str());
            if (key == "x_begin") xb = (int)val;
            else if (key == "x_end") xe = (int)val;
            else if (key == "y_begin") yb = (int)val;
            else if (key == "y_end") ye = (int)val;
            else if (key == "density") rho = val;
            else if (key == "velocity_u") u_vel = val;
            else if (key == "velocity_v") v_vel = val;
            else if (key == "pressure") p = val;
        }

        // Консервативные переменные
        double rhou = rho * u_vel;
        double rhov = rho * v_vel;
        double energy = p / (gamma - 1.0) + 0.5 * rho * (u_vel * u_vel + v_vel * v_vel);

        // Заполнение с учетом смещения
        for (int x = xb; x <= xe; ++x) {
            for (int y = yb; y <= ye; ++y) {
                int tx = x + offset_x;
                int ty = y + offset_y;
                if (tx < u.size() && ty < u[0].size()) {
                    u[tx][ty][0] = rho;
                    u[tx][ty][1] = rhou;
                    u[tx][ty][2] = rhov;
                    u[tx][ty][3] = energy;
                }
            }
        }
    }
}