#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <locale>
#include <algorithm>
#include <iomanip>
#include <stdexcept>
#include "Matrix.h"
#include "SimplexSolver.h"
#include "SimplexEngine.h"

int GLOBAL_TRIGGER = 0;



int main() {
    setlocale(LC_ALL, "Russian");
    try {
        SimplexEngine engine("data.txt");
        engine.printData(engine.mate, engine.operators, engine.function, engine.restrictions, engine.positive_coefs, engine.function_constant);
        SimplexSolver simplex(engine.function, engine.restrictions, engine.operators, engine.taskType, engine.mate, &engine, engine.positive_coefs, engine.function_constant, 0, 1);
        simplex.Solve();
        SimplexSolver simplexDouble(engine.function, engine.restrictions, engine.operators, engine.taskType, engine.mate, &engine, engine.positive_coefs, engine.function_constant, 1, 0);
        simplexDouble.Solve();
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }
    return 0;
}