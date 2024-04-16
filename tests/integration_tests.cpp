#include <fstream>

#include <gtest/gtest.h>

#include <Simulation.h>
#include <Util.h>
#include <string>

TEST(IntegrationTests, Simulation) {
    const SimulationConfig cfg = {
        .width = 10,
        .height = 10,
        .depth = 10,
        .endTime = 1.2,
        .stepSize = 0.3,
        .cellSize = 1.0,
        .outputPrefix = "integrationtest",
        .preconditionerType = PreconditionerType::DIC,
        .disableFileOutput = true,
    };

    SimulationOutput out = simulate(cfg);

    const auto velocityRaw = out.velocityField->getRawValues();

    // for (size_t i = 0; i < velocityRaw.getSize(); i++) {
    //     std::cout << velocityRaw[i] << "\n";
    // }

    std::string filename =
        std::string(TEST_RESOURCE_DIR) + "/velocity_snapshot.txt";
    std::ifstream infile(filename);

    size_t i = 0;
    ScalarT truth;
    while (infile >> truth) {
        ASSERT_LT(i, velocityRaw.getSize())
            << "Snapshot file has more entries than the simulation output";

        ASSERT_NEAR(velocityRaw[i], truth, 0.001)
            << "Simulation output and snapshot file differ at i = " +
                   std::to_string(i);

        i++;
    }

    ASSERT_EQ(i, velocityRaw.getSize())
        << "Simulation file has more entries than the snapshot file";
}
