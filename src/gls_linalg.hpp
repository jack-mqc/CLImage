// Copyright (c) 2021-2022 Glass Imaging Inc.
// Author: Fabio Riccardi <fabio@glass-imaging.com>
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef gls_linalg_h
#define gls_linalg_h

#include <iostream>
#include <stdexcept>

namespace gls {

template<int M, int N> struct Matrix;

// ---- Vector Type ----
template <int N>
struct Vector : public std::array<float, N> {
    Vector() { }

    Vector(const float(&il)[N]) {
        std::copy(il, il + N, this->begin());
    }

    Vector(const std::vector<float>& v) {
        assert(v.size() == N);
        std::copy(v.begin(), v.end(), this->begin());
    }

    Vector(std::initializer_list<float> list) {
        assert(list.size() == N);
        std::copy(list.begin(), list.end(), this->begin());
    }

    template<int P, int Q>
    requires (P * Q == N)
    Vector(const Matrix<P, Q>& m) {
        const auto ms = m.span();
        std::copy(ms.begin(), ms.end(), this->begin());
    }

    // Cast to a const float*
    operator const float*() const {
        return this->data();
    }
};

// Vector - Vector Multiplication (component-wise)
template <int N>
inline Vector<N> operator * (const Vector<N>& a, const Vector<N>& b) {
    auto ita = a.begin();
    auto itb = b.begin();
    Vector<N> result;
    std::for_each(result.begin(), result.end(), [&](float &r){ r = *ita++ * *itb++; });
    return result;
}

// Vector - Vector Division (component-wise)
template <int N>
inline Vector<N> operator / (const Vector<N>& a, const Vector<N>& b) {
    auto ita = a.begin();
    auto itb = b.begin();
    Vector<N> result;
    std::for_each(result.begin(), result.end(), [&](float &r){ r = *ita++ / *itb++; });
    return result;
}

// Vector - Scalar Addition
template <int N>
inline Vector<N> operator + (const Vector<N>& v, const float a) {
    auto itv = v.begin();
    Vector<N> result;
    std::for_each(result.begin(), result.end(), [&a, &itv](float &r){ r = *itv++ + a; });
    return result;
}

// Vector - Scalar Addition (commutative)
template <int N>
inline Vector<N> operator + (const float a, const Vector<N>& v) {
    return v + a;
}

// Vector - Scalar Subtraction
template <int N>
inline Vector<N> operator - (const Vector<N>& v, const float a) {
    auto itv = v.begin();
    Vector<N> result;
    std::for_each(result.begin(), result.end(), [&a, &itv](float &r){ r = *itv++ - a; });
    return result;
}

// Scalar - Vector Subtraction
template <int N>
inline Vector<N> operator - (const float a, const Vector<N>& v) {
    auto itv = v.begin();
    Vector<N> result;
    std::for_each(result.begin(), result.end(), [&a, &itv](float &r){ r = a - *itv++; });
    return result;
}

// Vector - Scalar Multiplication
template <int N>
inline Vector<N> operator * (const Vector<N>& v, const float a) {
    auto itv = v.begin();
    Vector<N> result;
    std::for_each(result.begin(), result.end(), [&a, &itv](float &r){ r = *itv++ * a; });
    return result;
}

// Scalar - Vector Multiplication (commutative)
template <int N>
inline Vector<N> operator * (const float a, const Vector<N>& v) {
    return v * a;
}

// Vector - Scalar Division
template <int N>
inline Vector<N> operator / (const Vector<N>& v, const float a) {
    auto itv = v.begin();
    Vector<N> result;
    std::for_each(result.begin(), result.end(), [&a, &itv](float &r){ r = *itv++ / a; });
    return result;
}

// Scalar - Vector Division
template <int N>
inline Vector<N> operator / (const float a, const Vector<N>& v) {
    auto itv = v.begin();
    Vector<N> result;
    std::for_each(result.begin(), result.end(), [&a, &itv](float &r){ r = a / *itv++; });
    return result;
}

// ---- Matrix Type ----

template <int N, int M>
struct Matrix : public std::array<Vector<M>, N> {
    Matrix() {}

    Matrix(const Vector<N * M>& v) {
        std::copy(v.begin(), v.end(), span().begin());
    }

    Matrix(const float(&il)[N * M]) {
        std::copy(il, il + (N * M), span().begin());
    }

    Matrix(const std::array<float, M>(&il)[N]) {
        // This is safe, il is just an array of arrays
        std::copy((float *) il, (float *) il + (N * M), span().begin());
    }

    Matrix(const std::vector<float>& v) {
        assert(v.size() == N * M);
        std::copy(v.begin(), v.end(), span().begin());
    }

    Matrix(std::initializer_list<float> list) {
        assert(list.size() == N * M);
        std::copy(list.begin(), list.end(), span().begin());
    }

    Matrix(std::initializer_list<std::array<float, M>> list) {
        assert(list.size() == N);
        int row = 0;
        for (const auto& v : list) {
            std::copy(v.begin(), v.end(), span(row++).begin());
        }
    }

    // Matrix Raw Data
    std::span<float> span() {
        return std::span(&(*this)[0][0], N * M);
    }

    const std::span<const float> span() const {
        return std::span(&(*this)[0][0], N * M);
    }

    // Matrix Row Raw Data
    std::span<float> span(int row) {
        return std::span(&(*this)[row][0], M);
    }

    const std::span<const float> span(int row) const {
        return std::span(&(*this)[row][0], M);
    }

    // Cast to a const float*
    operator const float*() const {
        return span().data();
    }
};

template <int N, int M>
std::span<float> span(Matrix<N, M>& m) {
    return std::span(&m[0][0], N * M);
}

template <int N, int M>
const std::span<const float> span(const Matrix<N, M>& m) {
    return std::span(&m[0][0], N * M);
}

// Matrix Transpose
template<int N, int M>
inline Matrix<N, M> transpose(const Matrix<M, N>& m) {
    Matrix<N, M> result;
    for (int j = 0; j < M; j++) {
        for (int i = 0; i < N; i++) {
            result[i][j] = m[j][i];
        }
    }
    return result;
}

// General Matrix Multiplication
template <int N, int K, int M>
inline Matrix<M, N> operator * (const Matrix<M, K>& a, const Matrix<K, N>& b) {
    Matrix<M, N> result;
    const auto bt = transpose(b);
    for (int j = 0; j < M; j++) {
        for (int i = 0; i < N; i++) {
            result[j][i] = 0;
            for (int k = 0; k < K; k++) {
                result[j][i] += a[j][k] * bt[i][k];
            }
        }
    }
    return result;
}

// Matrix - Vector Multiplication
template <int M, int N>
inline Vector<M> operator * (const Matrix<M, N>& a, const Vector<N>& b) {
    const auto result = a * Matrix<N, 1> { b };
    return Vector<M>(result);
}

// Vector - Matrix Multiplication
template <int M, int N>
inline Vector<N> operator * (const Vector<M>& a, const Matrix<M, N>& b) {
    const auto result = Matrix<1, N> { a } * b;
    return Vector<N>(result);
}

// (Square) Matrix Division (Multiplication with Inverse)
template <int N>
inline Matrix<N, N> operator / (const Matrix<N, N>& a, const Matrix<N, N>& b) {
    return a * inverse(b);
}

// Iterate over the elements of the input and output matrices applying a Matrix-Matrix function
template<int N, int M>
inline Matrix<N, M> apply(const Matrix<M, N>& a, const Matrix<M, N>& b, float (*f)(const float& a, const float& b)) {
    Matrix<N, M> result;
    auto ita = span(a).begin();
    auto itb = span(b).begin();
    for (auto& r : span(result)) {
        r = f(*ita++, *itb++);
    }
    return result;
}

// Iterate over the elements of the input and output matrices applying a Matrix-Scalar function
template<int N, int M>
inline Matrix<N, M> apply(const Matrix<M, N>& a, float b, float (*f)(const float& a, float b)) {
    Matrix<N, M> result;
    auto ita = span(a).begin();
    for (auto& r : span(result)) {
        r = f(*ita++, b);
    }
    return result;
}

// Matrix-Scalar Multiplication
template <int N, int M>
inline Matrix<N, M> operator * (const Matrix<N, M>& a, const float b) {
    return apply(a, b, [](const float& a, float b) {
        return a * b;
    });
}

// Matrix-Scalar Division
template <int N, int M>
inline Matrix<N, M> operator / (const Matrix<N, M>& a, const float b) {
    return apply(a, b, [](const float& a, float b) {
        return a / b;
    });
}

// Matrix-Matrix Addition
template <int N, int M>
inline Matrix<N, M> operator + (const Matrix<N, M>& a, const Matrix<N, M>& b) {
    return apply(a, b, [](const float& a, const float& b) {
        return a + b;
    });
}

// Matrix-Scalar Addition
template <int N, int M>
inline Matrix<N, M> operator + (const Matrix<N, M>& a, const float b) {
    return apply(a, b, [](const float& a, float b) {
        return a + b;
    });
}

// Matrix-Matrix Subtraction
template <int N, int M>
inline Matrix<N, M> operator - (const Matrix<N, M>& a, const Matrix<N, M>& b) {
    return apply(a, b, [](const float& a, const float& b) {
        return a - b;
    });
}

// Matrix-Scalar Subtraction
template <int N, int M>
inline Matrix<N, M> operator - (const Matrix<N, M>& a, const float b) {
    return apply(a, b, [](const float& a, float b) {
        return a - b;
    });
}

// --- Matrix Inverse Support ---

// Cofactor Matrix
// https://en.wikipedia.org/wiki/Minor_(linear_algebra)#Inverse_of_a_matrix
template <int N1, int N2 = N1 - 1>
inline Matrix<N2, N2> cofactor(const Matrix<N1, N1>& m, int p, int q) {
    assert(p < N1 && q < N1);

    Matrix<N2, N2> result;

    // Looping for each element of the matrix
    int i = 0, j = 0;
    for (int row = 0; row < N1; row++) {
        for (int col = 0; col < N1; col++) {
            //  Copying into temporary matrix only those element
            //  which are not in given row and column
            if (row != p && col != q) {
                result[i][j++] = m[row][col];

                // Row is filled, so increase row index and
                // reset col index
                if (j == N1 - 1) {
                    j = 0;
                    i++;
                }
            }
        }
    }
    return result;
}

// Matrix Determinant using Laplace's Cofactor Expansion
// https://en.wikipedia.org/wiki/Minor_(linear_algebra)#Cofactor_expansion_of_the_determinant
template <int N>
inline float determinant(const Matrix<N, N>& m) {
    assert(N > 1);

    float sign = 1;
    float result = 0;
    // Iterate for each element of first row
    for (int f = 0; f < N; f++) {
        result += sign * m[0][f] * determinant(cofactor(m, 0, f));
        // terms are to be added with alternate sign
        sign = -sign;
    }
    return result;
}

// Matrix Determinant, Special case for size 1x1
template <>
inline float determinant(const Matrix<1, 1>& m) {
    return m[0][0];
}

// Matrix Adjoint (Tanspose of the Cofactor Matrix)
// https://en.wikipedia.org/wiki/Adjugate_matrix
template <int N>
inline Matrix<N, N> adjoint(const Matrix<N, N>& m) {
    assert(N > 1);

    Matrix<N, N> adj;

    float sign = 1;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            // sign of adj[j][i] positive if sum of row
            // and column indexes is even.
            sign = ((i + j) % 2 == 0) ? 1 : -1;

            // Interchanging rows and columns to get the
            // transpose of the cofactor matrix
            float d = determinant(cofactor(m, i, j));
            adj[j][i] = d != 0 ? sign * d : 0;
        }
    }
    return adj;
}

// Matrix Adjoint - Special case for size 1x1
template <>
inline Matrix<1, 1> adjoint(const Matrix<1, 1>& m) {
    return { 1 };
}

// Inverse Matrix: inverse(m) = adj(m)/det(m)
// https://en.wikipedia.org/wiki/Minor_(linear_algebra)#Inverse_of_a_matrix
template <int N>
inline Matrix<N, N> inverse(const Matrix<N, N>& m) {
    float det = determinant(m);
    if (det == 0) {
        throw std::range_error("null determinant");
    }

    Matrix<N, N> inverse;
    const auto adj = adjoint(m);
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++) inverse[i][j] = adj[i][j] / det;

    return inverse;
}

// From DCRaw (https://www.dechifro.org/dcraw/)
template <int size>
gls::Matrix<size, 3> pseudoinverse(const gls::Matrix<size, 3>& in) {
    gls::Matrix<3,6> work;

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 6; j++) {
            work[i][j] = j == i + 3;
        }
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < size; k++) work[i][j] += in[k][i] * in[k][j];
        }
    }
    for (int i = 0; i < 3; i++) {
        float num = work[i][i];
        for (int j = 0; j < 6; j++) work[i][j] /= num;
        for (int k = 0; k < 3; k++) {
            if (k == i) continue;
            num = work[k][i];
            for (int j = 0; j < 6; j++) work[k][j] -= work[i][j] * num;
        }
    }
    gls::Matrix<size, 3> out;
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                out[i][j] += work[j][k + 3] * in[i][k];
            }
        }
    }
    return out;
}

// --- Utility Functions ---

template <int N>
std::ostream& operator<<(std::ostream& os, const Vector<N>& v) {
    for (int i = 0; i < N; i++) {
        os << v[i];
        if (i < N - 1) {
            os << ", ";
        }
    }
    return os;
}

template <int N, int M>
std::ostream& operator<<(std::ostream& os, const Matrix<N, M>& m) {
    for (int j = 0; j < N; j++) {
        for (int i = 0; i < M; i++) {
            os << m[j][i];
            if (j < N - 1 || i < M -1) {
                os << ", ";
            }
        }
        if (j < N-1) {
            os << std::endl;
        }
    }
    return os;
}

}  // namespace gls

#endif /* gls_linalg_h */
