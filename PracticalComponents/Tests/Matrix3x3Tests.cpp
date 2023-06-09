#include "pch.h"
#include "../PracticalComponents/PCM.h"

using namespace PC;

TEST(Matrix3x3, TestMatrixMultiplication) {
	Matrix3x3int a({ {1, 2, 3}, {4, 5, 6}, {7, 8, 9} }); //error here
	Matrix3x3int b({ {10, 11, 12}, {13, 14, 15}, {16, 17, 18} });
	Matrix3x3int result = a * b;

	Matrix3x3int expected({ {84, 90, 96}, {201, 216, 231}, {318, 342, 366} });

	EXPECT_EQ(result, expected);
}

TEST(Matrix3x3, TestMatrixAssignment) {
	Matrix3x3int a({ {1, 2, 3}, {4, 5, 6}, {7, 8, 9} });
	Matrix3x3int b = a;

	EXPECT_EQ(b, a);
}

TEST(Matrix3x3, TestMatrixTranspose) {
	Matrix3x3int a({ {1, 2, 3}, {4, 5, 6}, {7, 8, 9} });
	Matrix3x3int result = a.transpose();

	Matrix3x3int expected({ {1, 4, 7}, {2, 5, 8}, {3, 6, 9} });

	EXPECT_EQ(result, expected);
}

TEST(Matrix3x3, TestMatrixDeterminant) {
	Matrix3x3int a({ {1, 2, 3}, {4, 5, 6}, {7, 8, 9} });
	auto result = a.determinant();

	auto expected = 0;

	EXPECT_EQ(result, expected);
}

TEST(Matrix3x3, Matrix3x3FromTranslation) {
    Vector2<float> v{ 3.0f, 4.0f };
    Matrix3x3<float> mat = Matrix3x3<float>::Matrix3x3FromTranslation(v);

    // Validate expected translation matrix values
    EXPECT_EQ(mat.GetValue(0, 0), 1.0f);
    EXPECT_EQ(mat.GetValue(0, 2), v.x);
    EXPECT_EQ(mat.GetValue(1, 1), 1.0f);
    EXPECT_EQ(mat.GetValue(1, 2), v.y);
    EXPECT_EQ(mat.GetValue(2, 2), 1.0f);
}

TEST(Matrix3x3, Matrix3x3FromRotation) {
    float rotationAngle = 90.0f; // 90 degrees
    Matrix3x3<float> mat = Matrix3x3<float>::Matrix3x3FromRotation(rotationAngle);

    // Validate expected rotation matrix values with tolerance
    float tolerance = 1e-5f;
    EXPECT_NEAR(mat.GetValue(0, 0), cos(rotationAngle * M_PI / 180.0f), tolerance);
    EXPECT_NEAR(mat.GetValue(0, 1), -sin(rotationAngle * M_PI / 180.0f), tolerance);
    EXPECT_NEAR(mat.GetValue(1, 0), sin(rotationAngle * M_PI / 180.0f), tolerance);
    EXPECT_NEAR(mat.GetValue(1, 1), cos(rotationAngle * M_PI / 180.0f), tolerance);
    EXPECT_EQ(mat.GetValue(2, 2), 1.0f);
}

TEST(Matrix3x3, Matrix3x3FromScale) {
    Vector2<float> scaleFactors{ 2.0f, 3.0f };
    Matrix3x3<float> mat = Matrix3x3<float>::Matrix3x3FromScale(scaleFactors);

    // Validate expected scale matrix values
    EXPECT_EQ(mat.GetValue(0, 0), scaleFactors.x);
    EXPECT_EQ(mat.GetValue(1, 1), scaleFactors.y);
    EXPECT_EQ(mat.GetValue(2, 2), 1.0f);
}