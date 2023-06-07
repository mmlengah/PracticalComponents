#include "pch.h"
#include "../PracticalComponents/PCM.h"

using namespace PC;

TEST(Vector2, OperatorPlus) {
    // Initialize two vectors
    Vector2int v1(3, 4);
    Vector2int v2(1, 2);

    // Use the operator+
    Vector2int result = v1 + v2;

    // Assert the expected result
    EXPECT_EQ(result.x, 4);
    EXPECT_EQ(result.y, 6);
}