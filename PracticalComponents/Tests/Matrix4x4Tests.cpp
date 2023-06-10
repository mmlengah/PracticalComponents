#include "pch.h"
#include "../PracticalComponents/PCM.h"

using namespace PC;

TEST(Matrix4x4Test, EqualityTest) {
    Matrix4x4int m1 = { {1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11, 12}, {13, 14, 15, 16} };
    Matrix4x4int m2 = { {1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11, 12}, {13, 14, 15, 16} };

    ASSERT_TRUE(m1 == m2);
}

TEST(Matrix4x4Test, InequalityTest) {
    Matrix4x4int m1 = { {1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11, 12}, {13, 14, 15, 16} };
    Matrix4x4int m2 = { {1, 2, 3, 0}, {5, 6, 7, 8}, {9, 10, 11, 12}, {13, 14, 15, 16} };

    ASSERT_TRUE(m1 != m2);
}

TEST(Matrix4x4Test, MultiplicationTest) {
    Matrix4x4f m1 = { {1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11, 12}, {13, 14, 15, 16} };
    Matrix4x4f m2 = { {1, 1, 1, 1}, {1, 1, 1, 1}, {1, 1, 1, 1}, {1, 1, 1, 1} };

    Matrix4x4f m3 = m1 * m2;

    // The multiplication of m1 and m2 should result in each cell in m3 
    // being the sum of the corresponding row in m1 (because m2 is all 1's).
    for (int i = 0; i < 4; i++) {
        float rowSum = m1.GetValue(i, 0) + m1.GetValue(i, 1) + m1.GetValue(i, 2) + m1.GetValue(i, 3);
        for (int j = 0; j < 4; j++) {
            ASSERT_EQ(m3.GetValue(i, j), rowSum);
        }
    }
}
