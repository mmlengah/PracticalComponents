#include "pch.h"
#include "../PracticalComponents/PCM.h"

using namespace PC;

TEST(Vector3, OperatorPlus) {
    Vector3f a(1.0, 2.0, 3.0);
    Vector3f b(4.0, 5.0, 6.0);

    Vector3f result = a + b;

    EXPECT_FLOAT_EQ(result.x, 5.0);
    EXPECT_FLOAT_EQ(result.y, 7.0);
    EXPECT_FLOAT_EQ(result.z, 9.0);
}

TEST(Vector3, OperatorPlusEqual) {
    Vector3f a(1.0, 2.0, 3.0);
    Vector3f b(4.0, 5.0, 6.0);

    a += b;

    EXPECT_FLOAT_EQ(a.x, 5.0);
    EXPECT_FLOAT_EQ(a.y, 7.0);
    EXPECT_FLOAT_EQ(a.z, 9.0);
}

TEST(Vector3, OperatorPlusScalar) {
    Vector3f a(1.0, 2.0, 3.0);
    float scalar = 4.0f;

    Vector3f result = a + scalar;

    EXPECT_FLOAT_EQ(result.x, 5.0);
    EXPECT_FLOAT_EQ(result.y, 6.0);
    EXPECT_FLOAT_EQ(result.z, 7.0);
}

TEST(Vector3, OperatorPlusEqualScalar) {
    Vector3f a(1.0, 2.0, 3.0);
    float scalar = 4.0f;

    a += scalar;

    EXPECT_FLOAT_EQ(a.x, 5.0);
    EXPECT_FLOAT_EQ(a.y, 6.0);
    EXPECT_FLOAT_EQ(a.z, 7.0);
}

TEST(Vector3, OperatorMinus) {
    Vector3f a(4.0, 5.0, 6.0);
    Vector3f b(1.0, 2.0, 3.0);

    Vector3f result = a - b;

    EXPECT_FLOAT_EQ(result.x, 3.0);
    EXPECT_FLOAT_EQ(result.y, 3.0);
    EXPECT_FLOAT_EQ(result.z, 3.0);
}

TEST(Vector3, OperatorMinusEqual) {
    Vector3f a(4.0, 5.0, 6.0);
    Vector3f b(1.0, 2.0, 3.0);

    a -= b;

    EXPECT_FLOAT_EQ(a.x, 3.0);
    EXPECT_FLOAT_EQ(a.y, 3.0);
    EXPECT_FLOAT_EQ(a.z, 3.0);
}

TEST(Vector3, OperatorMinusScalar) {
    Vector3f a(4.0, 5.0, 6.0);
    float scalar = 1.0f;

    Vector3f result = a - scalar;

    EXPECT_FLOAT_EQ(result.x, 3.0);
    EXPECT_FLOAT_EQ(result.y, 4.0);
    EXPECT_FLOAT_EQ(result.z, 5.0);
}

TEST(Vector3, OperatorMinusEqualScalar) {
    Vector3f a(4.0, 5.0, 6.0);
    float scalar = 1.0f;

    a -= scalar;

    EXPECT_FLOAT_EQ(a.x, 3.0);
    EXPECT_FLOAT_EQ(a.y, 4.0);
    EXPECT_FLOAT_EQ(a.z, 5.0);
}

TEST(Vector3, OperatorMultiplyScalar) {
    Vector3f a(1.0, 2.0, 3.0);
    float scalar = 4.0f;

    Vector3f result = a * scalar;

    EXPECT_FLOAT_EQ(result.x, 4.0);
    EXPECT_FLOAT_EQ(result.y, 8.0);
    EXPECT_FLOAT_EQ(result.z, 12.0);
}

TEST(Vector3, OperatorMultiplyEqualScalar) {
    Vector3f a(1.0, 2.0, 3.0);
    float scalar = 4.0f;

    a *= scalar;

    EXPECT_FLOAT_EQ(a.x, 4.0);
    EXPECT_FLOAT_EQ(a.y, 8.0);
    EXPECT_FLOAT_EQ(a.z, 12.0);
}

TEST(Vector3, OperatorDivideScalar) {
    Vector3f a(4.0, 8.0, 12.0);
    float scalar = 4.0f;

    Vector3f result = a / scalar;

    EXPECT_FLOAT_EQ(result.x, 1.0);
    EXPECT_FLOAT_EQ(result.y, 2.0);
    EXPECT_FLOAT_EQ(result.z, 3.0);
}

TEST(Vector3, OperatorDivideEqualScalar) {
    Vector3f a(4.0, 8.0, 12.0);
    float scalar = 4.0f;

    a /= scalar;

    EXPECT_FLOAT_EQ(a.x, 1.0);
    EXPECT_FLOAT_EQ(a.y, 2.0);
    EXPECT_FLOAT_EQ(a.z, 3.0);
}

TEST(Vector3, OperatorEqual) {
    Vector3f a(1.0, 2.0, 3.0);
    Vector3f b(1.0, 2.0, 3.0);

    EXPECT_TRUE(a == b);
}

TEST(Vector3, OperatorNotEqual) {
    Vector3f a(1.0, 2.0, 3.0);
    Vector3f b(4.0, 5.0, 6.0);

    EXPECT_TRUE(a != b);
}

TEST(Vector3, OperatorLessThan) {
    Vector3f a(1.0, 2.0, 3.0); // magnitude sqrt(14)
    Vector3f b(4.0, 5.0, 6.0); // magnitude sqrt(77)

    EXPECT_TRUE(a < b);
}

TEST(Vector3, OperatorGreaterThan) {
    Vector3f a(4.0, 5.0, 6.0); // magnitude sqrt(77)
    Vector3f b(1.0, 2.0, 3.0); // magnitude sqrt(14)

    EXPECT_TRUE(a > b);
}

TEST(Vector3, OperatorLessThanEqual) {
    Vector3f a(1.0, 2.0, 3.0); // magnitude sqrt(14)
    Vector3f b(1.0, 2.0, 3.0); // magnitude sqrt(14)

    EXPECT_TRUE(a <= b);
}

TEST(Vector3, OperatorGreaterThanEqual) {
    Vector3f a(4.0, 5.0, 6.0); // magnitude sqrt(77)
    Vector3f b(4.0, 5.0, 6.0); // magnitude sqrt(77)

    EXPECT_TRUE(a >= b);
}

TEST(Vector3, DotProduct) {
    Vector3f a(1.0, 2.0, 3.0);
    Vector3f b(4.0, 5.0, 6.0);

    float dotProduct = a.dotProduct(b);

    EXPECT_FLOAT_EQ(dotProduct, 32.0f); // 1*4 + 2*5 + 3*6 = 32
}

TEST(Vector3, CrossProduct) {
    Vector3f a(1.0, 2.0, 3.0);
    Vector3f b(4.0, 5.0, 6.0);

    Vector3f crossProduct = a.crossProduct(b);

    EXPECT_FLOAT_EQ(crossProduct.x, -3.0);
    EXPECT_FLOAT_EQ(crossProduct.y, 6.0);
    EXPECT_FLOAT_EQ(crossProduct.z, -3.0);
}

TEST(Vector3, Magnitude) {
    Vector3f a(1.0, 2.0, 3.0);

    float magnitude = a.magnitude();

    EXPECT_FLOAT_EQ(magnitude, static_cast<float>(std::sqrt(14)));
}

TEST(Vector3, Normalise) {
    Vector3f a(1.0, 2.0, 3.0);

    Vector3f normalised = a.normalise();

    float magnitude = normalised.magnitude();
    EXPECT_FLOAT_EQ(magnitude, 1.0); // magnitude of normalised vector should be 1
}

TEST(Vector3, Angle) {
    Vector3f a(1.0, 0.0, 0.0);
    Vector3f b(0.0, 1.0, 0.0);

    float angle = a.angle(b);

    EXPECT_FLOAT_EQ(angle, 90.0); // Angle between orthogonal vectors should be 90 degrees
}

TEST(Vector3, Rotate) {
    Vector3f a(1.0, 0.0, 0.0);
    Vector3f axis(0.0, 0.0, 1.0);
    float angle = 90.0f;

    Vector3f rotated = a.rotate(axis, angle);

    EXPECT_NEAR(rotated.x, 0.0, 1e-4);
    EXPECT_NEAR(rotated.y, 1.0, 1e-4);
    EXPECT_NEAR(rotated.z, 0.0, 1e-4); // Assuming the right hand rule for rotation
}


