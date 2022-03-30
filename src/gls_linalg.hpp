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

template <int N, int M>
struct Matrix : std::array<std::array<float, M>, N> {
    std::span<float> span() {
        return std::span(&(*this)[0][0], N * M);
    }

    const std::span<const float> span() const {
        return std::span(&(*this)[0][0], N * M);
    }
};

typedef Matrix<3, 3> Mat;

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
            for (int k = 0; k < K; k++) {
                result[j][i] += a[j][k] * bt[i][k];
            }
        }
    }
    return result;
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
    auto ita = a.span().begin();
    auto itb = b.span().begin();
    for (auto& r : result.span()) {
        r = f(*ita++, *itb++);
    }
    return result;
}

// Iterate over the elements of the input and output matrices applying a Matrix-Scalar function
template<int N, int M>
inline Matrix<N, M> apply(const Matrix<M, N>& a, float b, float (*f)(const float& a, float b)) {
    Matrix<N, M> result;
    auto ita = a.span().begin();
    for (auto& r : result.span()) {
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

// --- Utility Functions ---

template <int N, int M>
inline void print(const Matrix<N, M>& m) {
    for (int j = 0; j < N; j++) {
        for (int i = 0; i < M; i++) {
            std::cout << m[j][i] << ", ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

template <int N, int M>
inline void print(const char* s, const Matrix<N, M>& m) {
    std::cout << s << ":" << std::endl;
    print(m);
}

}  // namespace gls

#endif /* gls_linalg_h */
