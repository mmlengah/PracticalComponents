#include "pch.h"
#include "../PracticalComponents/PCM.h"

using namespace PC;

TEST(Vector4, PlusOperatorOverloadWithVector) {
	Vector4int v1(1, 2, 3, 4);
	Vector4int v2(5, 6, 7, 8);

	Vector4int result = v1 + v2;

	EXPECT_EQ(result.x, 6);
	EXPECT_EQ(result.y, 8);
	EXPECT_EQ(result.z, 10);
	EXPECT_EQ(result.w, 12);
}

TEST(Vector4, PlusEqualsOperatorOverloadWithVector) {
	Vector4int v1(1, 2, 3, 4);
	Vector4int v2(5, 6, 7, 8);

	v1 += v2;

	EXPECT_EQ(v1.x, 6);
	EXPECT_EQ(v1.y, 8);
	EXPECT_EQ(v1.z, 10);
	EXPECT_EQ(v1.w, 12);
}

TEST(Vector4, PlusOperatorOverloadWithScalar) {
	Vector4int v1(1, 2, 3, 4);
	int scalar = 5;

	Vector4int result = v1 + scalar;

	EXPECT_EQ(result.x, 6);
	EXPECT_EQ(result.y, 7);
	EXPECT_EQ(result.z, 8);
	EXPECT_EQ(result.w, 9);
}

TEST(Vector4, PlusEqualsOperatorOverloadWithScalar) {
	Vector4int v1(1, 2, 3, 4);
	int scalar = 5;

	v1 += scalar;

	EXPECT_EQ(v1.x, 6);
	EXPECT_EQ(v1.y, 7);
	EXPECT_EQ(v1.z, 8);
	EXPECT_EQ(v1.w, 9);
}

TEST(Vector4, MinusOperatorOverloadWithVector) {
	Vector4int v1(5, 6, 7, 8);
	Vector4int v2(1, 2, 3, 4);

	Vector4int result = v1 - v2;

	EXPECT_EQ(result.x, 4);
	EXPECT_EQ(result.y, 4);
	EXPECT_EQ(result.z, 4);
	EXPECT_EQ(result.w, 4);
}

TEST(Vector4, MinusEqualsOperatorOverloadWithVector) {
	Vector4int v1(5, 6, 7, 8);
	Vector4int v2(1, 2, 3, 4);

	v1 -= v2;

	EXPECT_EQ(v1.x, 4);
	EXPECT_EQ(v1.y, 4);
	EXPECT_EQ(v1.z, 4);
	EXPECT_EQ(v1.w, 4);
}

TEST(Vector4, MinusOperatorOverloadWithScalar) {
	Vector4int v1(5, 6, 7, 8);
	int scalar = 1;

	Vector4int result = v1 - scalar;

	EXPECT_EQ(result.x, 4);
	EXPECT_EQ(result.y, 5);
	EXPECT_EQ(result.z, 6);
	EXPECT_EQ(result.w, 7);
}

TEST(Vector4, MinusEqualsOperatorOverloadWithScalar) {
	Vector4int v1(5, 6, 7, 8);
	int scalar = 1;

	v1 -= scalar;

	EXPECT_EQ(v1.x, 4);
	EXPECT_EQ(v1.y, 5);
	EXPECT_EQ(v1.z, 6);
	EXPECT_EQ(v1.w, 7);
}

TEST(Vector4, MultiplyOperatorOverloadWithScalar) {
	Vector4int v1(1, 2, 3, 4);
	int scalar = 2;

	Vector4int result = v1 * scalar;

	EXPECT_EQ(result.x, 2);
	EXPECT_EQ(result.y, 4);
	EXPECT_EQ(result.z, 6);
	EXPECT_EQ(result.w, 8);
}

TEST(Vector4, MultiplyEqualsOperatorOverloadWithScalar) {
	Vector4int v1(1, 2, 3, 4);
	int scalar = 2;

	v1 *= scalar;

	EXPECT_EQ(v1.x, 2);
	EXPECT_EQ(v1.y, 4);
	EXPECT_EQ(v1.z, 6);
	EXPECT_EQ(v1.w, 8);
}

TEST(Vector4, DivideOperatorOverloadWithScalar) {
	Vector4int v1(2, 4, 6, 8);
	int scalar = 2;

	Vector4int result = v1 / scalar;

	EXPECT_EQ(result.x, 1);
	EXPECT_EQ(result.y, 2);
	EXPECT_EQ(result.z, 3);
	EXPECT_EQ(result.w, 4);
}

TEST(Vector4, DivideEqualsOperatorOverloadWithScalar) {
	Vector4int v1(2, 4, 6, 8);
	int scalar = 2;

	v1 /= scalar;

	EXPECT_EQ(v1.x, 1);
	EXPECT_EQ(v1.y, 2);
	EXPECT_EQ(v1.z, 3);
	EXPECT_EQ(v1.w, 4);
}

TEST(Vector4, DivisionByZero) {
	Vector4int v1(1, 2, 3, 4);
	int scalar = 0;

	ASSERT_DEATH(v1 / scalar, ".*");  // Assuming that the program exits when a division by zero occurs
	ASSERT_DEATH(v1 /= scalar, ".*"); // We're using the ASSERT_DEATH function to check if the program terminates as expected
}

TEST(Vector4, EqualityOperator) {
	Vector4int v1(1, 2, 3, 4);
	Vector4int v2(1, 2, 3, 4);

	EXPECT_TRUE(v1 == v2);
}

TEST(Vector4, InequalityOperator) {
	Vector4int v1(1, 2, 3, 4);
	Vector4int v2(5, 6, 7, 8);

	EXPECT_TRUE(v1 != v2);
}

TEST(Vector4, DotProduct) {
	Vector4int v1(1, 2, 3, 4);
	Vector4int v2(5, 6, 7, 8);

	EXPECT_EQ(v1.dotProduct(v2), 70); // 1*5 + 2*6 + 3*7 + 4*8 = 70
}

TEST(Vector4, Magnitude) {
	Vector4int v1(1, 2, 3, 4);

	EXPECT_EQ(v1.magnitude(), static_cast<int>(std::sqrt(30))); // sqrt(1^2 + 2^2 + 3^2 + 4^2) = sqrt(30)
}

TEST(Vector4, Normalise) {
	Vector4int v1(1, 2, 3, 4);
	Vector4int result = v1.normalise();
	int mag = v1.magnitude();

	EXPECT_EQ(result.x, static_cast<int>(1 / mag));
	EXPECT_EQ(result.y, static_cast<int>(2 / mag));
	EXPECT_EQ(result.z, static_cast<int>(3 / mag));
	EXPECT_EQ(result.w, static_cast<int>(4 / mag));
}

