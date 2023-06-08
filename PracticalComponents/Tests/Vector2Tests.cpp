#include "pch.h"
#include "../PracticalComponents/PCM.h"

using namespace PC;

// Test the addition operator +
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

// Test the addition operator +=
TEST(Vector2, OperatorPlusEqualsTest) {
    Vector2<int> vec1(1, 2);
    Vector2<int> vec2(3, 4);
    vec1 += vec2;
    EXPECT_EQ(vec1.x, 4);
    EXPECT_EQ(vec1.y, 6);
}

// Test the scalar addition operator +
TEST(Vector2, OperatorPlusScalarTest) {
    Vector2<int> vec(1, 2);
    int scalar = 3;
    Vector2<int> result = vec + scalar;
    EXPECT_EQ(result.x, 4);
    EXPECT_EQ(result.y, 5);
}

// Test the scalar addition operator +=
TEST(Vector2, OperatorPlusEqualsScalarTest) {
    Vector2<int> vec(1, 2);
    int scalar = 3;
    vec += scalar;
    EXPECT_EQ(vec.x, 4);
    EXPECT_EQ(vec.y, 5);
}

// Test the subtraction operator -
TEST(Vector2, OperatorMinusTest) {
    Vector2<int> vec1(1, 2);
    Vector2<int> vec2(3, 4);
    Vector2<int> result = vec1 - vec2;
    EXPECT_EQ(result.x, -2);
    EXPECT_EQ(result.y, -2);
}

// Test the subtraction operator -=
TEST(Vector2, OperatorMinusEqualsTest) {
    Vector2<int> vec1(1, 2);
    Vector2<int> vec2(3, 4);
    vec1 -= vec2;
    EXPECT_EQ(vec1.x, -2);
    EXPECT_EQ(vec1.y, -2);
}

// Test the scalar multiplication operator *
TEST(Vector2, OperatorMultiplyScalarTest) {
    Vector2<int> vec(2, 3);
    int scalar = 3;
    Vector2<int> result = vec * scalar;
    EXPECT_EQ(result.x, 6);
    EXPECT_EQ(result.y, 9);
}

// Test the scalar multiplication operator *=
TEST(Vector2, OperatorMultiplyEqualsScalarTest) {
    Vector2<int> vec(2, 3);
    int scalar = 3;
    vec *= scalar;
    EXPECT_EQ(vec.x, 6);
    EXPECT_EQ(vec.y, 9);
}

// Test the scalar division operator /
TEST(Vector2, OperatorDivideScalarTest) {
    Vector2<int> vec(6, 9);
    int scalar = 3;
    Vector2<int> result = vec / scalar;
    EXPECT_EQ(result.x, 2);
    EXPECT_EQ(result.y, 3);
}

// Test the scalar division operator /=
TEST(Vector2, OperatorDivideEqualsScalarTest) {
    Vector2<int> vec(6, 9);
    int scalar = 3;
    vec /= scalar;
    EXPECT_EQ(vec.x, 2);
    EXPECT_EQ(vec.y, 3);
}

// Test the equality operator ==
TEST(Vector2, OperatorEqualsTest) {
    Vector2<int> vec1(2, 3);
    Vector2<int> vec2(2, 3);
    EXPECT_TRUE(vec1 == vec2);
}

// Test the inequality operator !=
TEST(Vector2, OperatorNotEqualsTest) {
    Vector2<int> vec1(2, 3);
    Vector2<int> vec2(3, 2);
    EXPECT_TRUE(vec1 != vec2);
}

// Test the unary minus operator -
TEST(Vector2, UnaryMinusTest) {
    Vector2<int> vec(2, 3);
    Vector2<int> result = -vec;
    EXPECT_EQ(result.x, -2);
    EXPECT_EQ(result.y, -3);
}

// Test the less than operator <
TEST(Vector2, OperatorLessThanTest) {
    Vector2<int> vec1(2, 2);
    Vector2<int> vec2(3, 3); // The magnitude of vec2 will be greater than the magnitude of vec1.
    EXPECT_TRUE(vec1 < vec2);
}

// Test the greater than operator >
TEST(Vector2, OperatorGreaterThanTest) {
    Vector2<int> vec1(3, 3);
    Vector2<int> vec2(2, 2); // The magnitude of vec1 will be greater than the magnitude of vec2.
    EXPECT_TRUE(vec1 > vec2);
}

// Test the less than or equal operator <=
TEST(Vector2, OperatorLessThanOrEqualTest) {
    Vector2<int> vec1(2, 2);
    Vector2<int> vec2(2, 2); // The magnitudes of vec1 and vec2 are equal.
    Vector2<int> vec3(3, 3); // The magnitude of vec3 is greater than the magnitude of vec1.
    EXPECT_TRUE(vec1 <= vec2);
    EXPECT_TRUE(vec1 <= vec3);
}

// Test the greater than or equal operator >=
TEST(Vector2, OperatorGreaterThanOrEqualTest) {
    Vector2<int> vec1(3, 3);
    Vector2<int> vec2(2, 2); // The magnitude of vec1 is greater than the magnitude of vec2.
    Vector2<int> vec3(3, 3); // The magnitudes of vec1 and vec3 are equal.
    EXPECT_TRUE(vec1 >= vec2);
    EXPECT_TRUE(vec1 >= vec3);
}

// Test dotProduct method
TEST(Vector2, DotProductTest) {
    Vector2<int> vec1(2, 3);
    Vector2<int> vec2(3, 4);
    int dotProd = vec1.dotProduct(vec2);
    EXPECT_EQ(dotProd, 18); // 2*3 + 3*4
}

// Test magnitude method
TEST(Vector2, MagnitudeTest) {
    Vector2<int> vec(3, 4);
    int mag = vec.magnitude();
    EXPECT_EQ(mag, 5); // sqrt(3*3 + 4*4)
}

// Test normalize method
TEST(Vector2Fl, NormalizeTest) {
    Vector2<float> vec(3.0f, 4.0f);
    Vector2<float> norm = vec.normalize();
    EXPECT_NEAR(norm.x, 0.6f, 0.0001f);
    EXPECT_NEAR(norm.y, 0.8f, 0.0001f);
}

// Test angle method
TEST(Vector2Fl, AngleTest) {
    Vector2<float> vec1(1.0f, 0.0f);
    Vector2<float> vec2(0.0f, 1.0f);
    float angle = vec1.angle(vec2);
    EXPECT_NEAR(angle, 90.0f, 0.0001f);
}

// Test rotate method
TEST(Vector2Fl, RotateTest) {
    Vector2<float> vec(1.0f, 0.0f);
    Vector2<float> rotated = vec.rotate(90.0f);
    EXPECT_NEAR(rotated.x, 0.0f, 0.0001f);
    EXPECT_NEAR(rotated.y, 1.0f, 0.0001f);
}
