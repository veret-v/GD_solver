#ifndef SCRIPT
#define SCRIPT

#include <vector>
#include <string>

// Добавляем параметры для учета фиктивных ячеек и физических констант
void parseAndInitialize(std::vector<std::vector<std::vector<double>>>& u, 
                        const std::string& filename, 
                        int offset_x, 
                        int offset_y, 
                        double gamma);

#endif