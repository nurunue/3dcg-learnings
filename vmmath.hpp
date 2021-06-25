//! @file vmmath.hpp
//! @brief 3DCGに使用する基本的な数学
//! @author suzulang.com
//! @date 令和三年六月二十五日

#pragma once

#include <cmath>
namespace szl {
	namespace maths {

		//! @brief 入力が浮動小数点の時はそのまま返し、それ以外の時はfloatを返す
		template<typename T, bool isint>struct IsRealT;
		template<typename T>struct IsRealT<T, true> { using type = T; };
		template<typename T>struct IsRealT<T, false> { using type = float; };

		//! @brief 指定した型が整数ならfloatを返すテンプレート
		template<typename T>
		struct RealIfInt {
			using type = typename IsRealT<T, std::is_floating_point<T>::value>::type;
		};

		//! @brief 誤差の定義
		template<typename T>
		struct real_error { static constexpr float value = (float)1e-7; };
		template<>
		struct real_error<double> { static constexpr double value = 1e-15; };


	}
	namespace maths {

		//! @brief 誤差の閾値未満かをチェック
		//! @param [in] value 検査する値
		//! @retval true 値の絶対値は閾値より小さい（valueは0にとても近い）
		//! @retval false 値の絶対値は閾値より大きい（valueは0ではない）
		template<typename T>
		inline bool is_error(const T value) {
			return
				abs(value) < real_error<T>::value;
		}

		//! @brief 与えたラジアンを度に変換する
		//! @param [in] radian ラジアンの値
		//! @return 度の値
		template<typename ScalarT, typename ReturnT = typename maths::RealIfInt<ScalarT>::type>
		inline auto to_degree(const ScalarT radian) {
			return static_cast<ReturnT>(radian * 180.0 / 3.1415926535897932384626);
		}

		//! @brief 与えた度をラジアンに変換する
		//! @param [in] degree 度の値
		//! @return ラジアンの値
		template<typename ScalarT, typename ReturnT = typename maths::RealIfInt<ScalarT>::type>
		inline auto to_radian(const ScalarT degree) {
			return static_cast<ReturnT>(degree * 3.1415926535897932384626 / 180.0);
		}


		//! @brief 値を指定した範囲に収める
		//! @param [in] value
		//! @param [in] low valueが取り得る最小値
		//! @param [in] high valueが取り得る最大値
		//! @return low <= value <= high が保証されたvalue
		template<typename VALUE, typename LOW, typename HIGH>
		inline VALUE clamp(const VALUE value, const LOW low, const HIGH high) {
			if (value < low)
				return static_cast<VALUE>(low);
			if (value > high)
				return static_cast<VALUE>(high);
			return value;
		}

		//! @brief cotangent
		template<typename ScalarT, typename ReturnT = typename maths::RealIfInt<ScalarT>::type>
		inline ReturnT cot(const ScalarT theta) {
			return (ReturnT)(1.0 / tan(theta));
		}

	}

	//ベクトル
	namespace mathv {

		//! @brief 三次元ベクトルの長さを求める
		//! @param [in] v 三次元ベクトル
		//! @return ベクトルの長さ
		template<typename VectorT3>
		inline auto length3(const VectorT3& v) {
			return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
		}

		//! @brief ベクトルの足し算
		//! @param [in] vd 結果の格納先
		//! @param [in] va 三次元ベクトル
		//! @param [in] vb 三次元ベクトル
		//! @return vd
		template<typename VectorTD3, typename VectorTA3, typename VectorTB3>
		inline VectorTD3& add3(VectorTD3& vd, const VectorTA3& va, const VectorTB3& vb) {

			using TDScalarT = typename std::remove_reference<decltype(vd[0])>::type;

			const auto x1 = va[0];
			const auto y1 = va[1];
			const auto z1 = va[2];
			const auto x2 = vb[0];
			const auto y2 = vb[1];
			const auto z2 = vb[2];

			vd[0] = static_cast<TDScalarT>(va[0] + vb[0]);
			vd[1] = static_cast<TDScalarT>(va[1] + vb[1]);
			vd[2] = static_cast<TDScalarT>(va[2] + vb[2]);

			return vd;

		}
		//ベクトルの引き算
		//! @param [out] vd 結果の格納先
		//! @param [in] va 三次元ベクトル
		//! @param [in] vb 三次元ベクトル
		//! @return vd
		template<typename VectorTD3, typename VectorTA3, typename VectorTB3>
		inline VectorTD3& sub3(VectorTD3& vd, const VectorTA3& va, const VectorTB3& vb) {

			using TDScalarT = typename std::remove_reference<decltype(vd[0])>::type;

			vd[0] = static_cast<TDScalarT>(va[0] - vb[0]);
			vd[1] = static_cast<TDScalarT>(va[1] - vb[1]);
			vd[2] = static_cast<TDScalarT>(va[2] - vb[2]);

			return vd;

		}

		//! @brief 三次元ベクトルの内積を求める
		//! @param [in] va 三次元ベクトル
		//! @param [in] vb 三次元ベクトル
		//! @return 内積
		template<typename VectorTA3, typename VectorTB3>
		inline auto inner3(const VectorTA3& va, const VectorTB3& vb) {
			return ((va[0])*(vb[0]) + (va[1])*(vb[1]) + (va[2])*(vb[2]));
		}

		//! @brief 三次元ベクトルの外積を求める
		//! @param [out] vd 結果の格納先
		//! @param [in] va 三次元ベクトル
		//! @param [in] vb 三次元ベクトル
		//! @return vd
		template<typename VectorTD3, typename VectorTA3, typename VectorTB3>
		inline VectorTD3& outer3(VectorTD3& vd, const VectorTA3& va, const VectorTB3& vb) {

			using TDScalarT = typename std::remove_reference<decltype(vd[0])>::type;

			const auto x1 = va[0];
			const auto y1 = va[1];
			const auto z1 = va[2];
			const auto x2 = vb[0];
			const auto y2 = vb[1];
			const auto z2 = vb[2];

			vd[0] = static_cast<TDScalarT>(y1 * z2 - z1 * y2);
			vd[1] = static_cast<TDScalarT>(z1 * x2 - x1 * z2);
			vd[2] = static_cast<TDScalarT>(x1 * y2 - y1 * x2);

			return vd;

		}

		//! @brief 三次元ベクトルを正規化する
		//! @param [in] v 三次元ベクトル
		//! @retval true 計算成功
		//! @retval false ベクトルの長さが短すぎて失敗した
		template<typename VectorT3>
		inline bool normalize3(VectorT3& v)
		{
			auto& x = v[0];
			auto& y = v[1];
			auto& z = v[2];

			auto len = std::sqrt(x * x + y * y + z * z);

			if (maths::is_error(len) == true)
				return false;

			len = static_cast<decltype(len)>(1.0) / len;
			x *= len;
			y *= len;
			z *= len;

			return true;
		}

		//! @brief 三次元ベクトルの角度を計算する
		//! @param [in] va 三次元ベクトル
		//! @param [in] vb 三次元ベクトル
		//! @return ラジアンの角度
		template<typename VectorTA3, typename VectorTB3>
		inline auto angle3(const VectorTA3& va, const VectorTB3 vb) {
			auto ac = inner3(va, vb) / (length3(va) * length3(vb));
			using scalarT = decltype(ac);
			constexpr auto low = static_cast<scalarT>(-1.0) + maths::real_error<scalarT>::value;
			constexpr auto high = static_cast<scalarT>(1.0) - maths::real_error<scalarT>::value;
			return std::acos(clamp(ac, low, high));
		}


		//! @brief 直線の始点・終点から、始点→終点のベクトルを求める
		//! @param [out] vd 結果の格納先
		//! @param [in] from 直線の始点
		//! @param [in] to 直線の終点
		//! @return vd
		template<typename VectorTD3, typename VectorFrom3, typename VectorTo3>
		inline VectorTD3& vector_from_line3(VectorTD3& vd, const VectorFrom3& from, const VectorTo3& to) {
			vd[0] = to[0] - from[0];
			vd[1] = to[1] - from[1];
			vd[2] = to[2] - from[2];
			return vd;
		}
		//! @brief 二つの座標の距離
		//! @param [in] va 三次元座標１
		//! @param [in] vb 三次元座標２
		//! return 距離
		template<typename VectorTA3, typename VectorTB3>
		inline auto distance3(const VectorTA3& va, const VectorTB3& vb) {
			return sqrt(
				pow((vb[0] - va[0]), 2) +
				pow((vb[1] - va[1]), 2) +
				pow((vb[2] - va[2]), 2));
		}

		//! @brief ベクトルにスカラーをかける
		//! @param [in,out] vd 元となるベクトル及び計算結果
		//! @param [in] scalar 掛けるスカラー値
		//! @return vd
		template<typename VectorT3, typename ScalarT>
		inline VectorT3& multi3(VectorT3& vd, const ScalarT scalar) {
			vd[0] *= scalar;
			vd[1] *= scalar;
			vd[2] *= scalar;
			return vd;
		}

		//! @brief 配列にNANが入っているか確認
		//! @param [in] v 評価する配列
		//! @param [in] count 行列の要素数
		//! @retval true 配列の中にnanがある
		//! @retval false 配列の中にnanはない
		template<typename ScalarT>
		bool chek_array_nan(const ScalarT& v, const size_t count) {
			for (size_t i = 0; i < count; i++) {
				if (isnan(v[i]) == true)
					return true;
			}
			return false;
		}

	}


	//行列（補助関数）
	namespace mathm {
		////////////////////////////////////////////////
		// 要素番号計算
		////////////////////////////////////////////////

		//! @brief 与えられた番号が正方行列の対角線の要素かどうかを判定する
		//! @param [in] index 一次元配列の要素番号
		//! @param [in] Row 行列の一行の長さ
		//! @retval indexが一行Rowの行列上でi==j
		//! @retval i != j
		inline bool is_identity_matrix_diagonal(const int index, const int Row) {
			int row = index / Row;
			int col = index % Row;
			return row == col;
		}

		//! @brief IxJ行列の i,j の要素にアクセスする
		//! @param [in] i 行番号
		//! @param [in] j 列番号
		//! @param [in] I 行数
		//! @return i,jの一次元配列におけるインデクス
		inline int matijI(const int i, const int j, const int I) {
			return i + I * j;
		}
	}

	// 行列演算用関数群
	namespace mathm {

		////////////////////////////////////////////////
		// 零行列
		////////////////////////////////////////////////
		//! @brief 行列の要素を全て０にする
		//! @param [in] _mat operator[]で一次元配列としてアクセス可能な行列
		//! @param [in] size size一次元配列としての要素数
		//! @return なし
		template<typename Matrix>
		void zero_matrix(Matrix& _mat, const int size) {
			using TT = typename std::remove_reference<decltype(_mat[0])>::type;
			for (size_t i = 0; i < size; i++) {
				_mat[i] = TT(0);
			}
		}

		//! @brief 行列の要素を全て０にする(4x4行列用)
		//! @return なし
		template<typename Matrix>
		inline void zero_matrix4(Matrix& mat44) {
			zero_matrix(mat44, 4 * 4);
		}
		//! @brief 行列の要素を全て０にする(3x3行列用)
		//! @return なし
		template<typename Matrix>
		inline void zero_matrix3(Matrix& mat44) {
			zero_matrix(mat44, 3 * 3);
		}

		////////////////////////////////////////////////
		// 単位行列
		////////////////////////////////////////////////

		//! @brief 正方行列を単位行列にする
		//! @param [in,out] _mat 正方行列
		//! @param [in] Row 一行の要素数
		//! @return _mat
		template<typename Matrix>
		Matrix& load_identity(Matrix& _mat, const int Row) {
			int size = Row * Row;
			for (int i = 0; i < size; i++) {
				if (is_identity_matrix_diagonal(i, Row)) {
					_mat[i] = 1;
				}
				else {
					_mat[i] = 0;
				}
			}
			return _mat;
		}

		//! @brief 正方行列を単位行列にする(4x4)
		//! @param [in,out] _mat 正方行列
		//! @return _mat
		template<typename Matrix>
		Matrix& load_identity4(Matrix& _mat) {
			return load_identity(_mat, 4);
		}

		//! @brief 正方行列を単位行列にする(3x3)
		//! @param [in,out] _mat 正方行列
		//! @return _mat
		template<typename Matrix>
		Matrix& load_identity3(Matrix& _mat) {
			return load_identity(_mat, 3);
		}


		////////////////////////////////////////////////
		// 逆行列
		////////////////////////////////////////////////

		//! @brief 4x4行列の行列式を取得
		//! @param [in] dm16 4x4行列
		//! @return 行列式
		template<typename Matrix>
		inline auto determinant44(const Matrix& dm16) {

			using scalar_t = typename std::remove_reference<decltype(dm16[0])>::type;

			const scalar_t a11 = dm16[0];
			const scalar_t a12 = dm16[1];
			const scalar_t a13 = dm16[2];
			const scalar_t a14 = dm16[3];

			const scalar_t a21 = dm16[4];
			const scalar_t a22 = dm16[5];
			const scalar_t a23 = dm16[6];
			const scalar_t a24 = dm16[7];

			const scalar_t a31 = dm16[8];
			const scalar_t a32 = dm16[9];
			const scalar_t a33 = dm16[10];
			const scalar_t a34 = dm16[11];

			const scalar_t a41 = dm16[12];
			const scalar_t a42 = dm16[13];
			const scalar_t a43 = dm16[14];
			const scalar_t a44 = dm16[15];


			return a11 * a22* a33* a44 + a11 * a23 * a34 * a42 + a11 * a24 * a32 * a43
				+ a12 * a21 * a34 * a43 + a12 * a23 * a31 * a44 + a12 * a24 * a33 * a41
				+ a13 * a21 * a32 * a44 + a13 * a22 * a34 * a41 + a13 * a24 * a31 * a42
				+ a14 * a21 * a33 * a42 + a14 * a22 * a31 * a43 + a14 * a23 * a32 * a41
				- a11 * a22 * a34 * a43 - a11 * a23 * a32 * a44 - a11 * a24 * a33 * a42
				- a12 * a21 * a33 * a44 - a12 * a23 * a34 * a41 - a12 * a24 * a31 * a43
				- a13 * a21 * a34 * a42 - a13 * a22 * a31 * a44 - a13 * a24 * a32 * a41
				- a14 * a21 * a32 * a43 - a14 * a22 * a33 * a41 - a14 * a23 * a31 * a42;
		}





		//! @brief 4x4行列の逆行列を取得
		//! @param [out] dst_dm16 結果を格納する一次元配列
		//! @param [in] src_dm16 元の行列
		//! @retval true 成功した
		//! @retval false 行列式の大きさが小さすぎて失敗した
		template<typename MatrixD, typename MatrixS>
		inline bool inverse44(MatrixD& dst_dm16, const MatrixS& src_dm16) {

			using d_scalar_t = typename std::remove_reference<decltype(dst_dm16[0])>::type;
			using s_scalar_t = typename std::remove_reference<decltype(src_dm16[0])>::type;

			s_scalar_t a11 = src_dm16[0];
			s_scalar_t a12 = src_dm16[1];
			s_scalar_t a13 = src_dm16[2];
			s_scalar_t a14 = src_dm16[3];

			s_scalar_t a21 = src_dm16[4];
			s_scalar_t a22 = src_dm16[5];
			s_scalar_t a23 = src_dm16[6];
			s_scalar_t a24 = src_dm16[7];

			s_scalar_t a31 = src_dm16[8];
			s_scalar_t a32 = src_dm16[9];
			s_scalar_t a33 = src_dm16[10];
			s_scalar_t a34 = src_dm16[11];

			s_scalar_t a41 = src_dm16[12];
			s_scalar_t a42 = src_dm16[13];
			s_scalar_t a43 = src_dm16[14];
			s_scalar_t a44 = src_dm16[15];

			d_scalar_t& b11 = dst_dm16[0];
			d_scalar_t& b12 = dst_dm16[1];
			d_scalar_t& b13 = dst_dm16[2];
			d_scalar_t& b14 = dst_dm16[3];
			d_scalar_t& b21 = dst_dm16[4];
			d_scalar_t& b22 = dst_dm16[5];

			d_scalar_t& b23 = dst_dm16[6];
			d_scalar_t& b24 = dst_dm16[7];
			d_scalar_t& b31 = dst_dm16[8];
			d_scalar_t& b32 = dst_dm16[9];
			d_scalar_t& b33 = dst_dm16[10];
			d_scalar_t& b34 = dst_dm16[11];

			d_scalar_t& b41 = dst_dm16[12];
			d_scalar_t& b42 = dst_dm16[13];
			d_scalar_t& b43 = dst_dm16[14];
			d_scalar_t& b44 = dst_dm16[15];

			b11 = static_cast<d_scalar_t>(a22 * a33 * a44 + a23 * a34 * a42 + a24 * a32 * a43 - a22 * a34 * a43 - a23 * a32 * a44 - a24 * a33 * a42);
			b12 = static_cast<d_scalar_t>(a12 * a34 * a43 + a13 * a32 * a44 + a14 * a33 * a42 - a12 * a33 * a44 - a13 * a34 * a42 - a14 * a32 * a43);
			b13 = static_cast<d_scalar_t>(a12 * a23 * a44 + a13 * a24 * a42 + a14 * a22 * a43 - a12 * a24 * a43 - a13 * a22 * a44 - a14 * a23 * a42);
			b14 = static_cast<d_scalar_t>(a12 * a24 * a33 + a13 * a22 * a34 + a14 * a23 * a32 - a12 * a23 * a34 - a13 * a24 * a32 - a14 * a22 * a33);

			b21 = static_cast<d_scalar_t>(a21 * a34 * a43 + a23 * a31 * a44 + a24 * a33 * a41 - a21 * a33 * a44 - a23 * a34 * a41 - a24 * a31 * a43);
			b22 = static_cast<d_scalar_t>(a11 * a33 * a44 + a13 * a34 * a41 + a14 * a31 * a43 - a11 * a34 * a43 - a13 * a31 * a44 - a14 * a33 * a41);
			b23 = static_cast<d_scalar_t>(a11 * a24 * a43 + a13 * a21 * a44 + a14 * a23 * a41 - a11 * a23 * a44 - a13 * a24 * a41 - a14 * a21 * a43);
			b24 = static_cast<d_scalar_t>(a11 * a23 * a34 + a13 * a24 * a31 + a14 * a21 * a33 - a11 * a24 * a33 - a13 * a21 * a34 - a14 * a23 * a31);
			b31 = static_cast<d_scalar_t>(a21 * a32 * a44 + a22 * a34 * a41 + a24 * a31 * a42 - a21 * a34 * a42 - a22 * a31 * a44 - a24 * a32 * a41);
			b32 = static_cast<d_scalar_t>(a11 * a34 * a42 + a12 * a31 * a44 + a14 * a32 * a41 - a11 * a32 * a44 - a12 * a34 * a41 - a14 * a31 * a42);
			b33 = static_cast<d_scalar_t>(a11 * a22 * a44 + a12 * a24 * a41 + a14 * a21 * a42 - a11 * a24 * a42 - a12 * a21 * a44 - a14 * a22 * a41);
			b34 = static_cast<d_scalar_t>(a11 * a24 * a32 + a12 * a21 * a34 + a14 * a22 * a31 - a11 * a22 * a34 - a12 * a24 * a31 - a14 * a21 * a32);
			b41 = static_cast<d_scalar_t>(a21 * a33 * a42 + a22 * a31 * a43 + a23 * a32 * a41 - a21 * a32 * a43 - a22 * a33 * a41 - a23 * a31 * a42);
			b42 = static_cast<d_scalar_t>(a11 * a32 * a43 + a12 * a33 * a41 + a13 * a31 * a42 - a11 * a33 * a42 - a12 * a31 * a43 - a13 * a32 * a41);
			b43 = static_cast<d_scalar_t>(a11 * a23 * a42 + a12 * a21 * a43 + a13 * a22 * a41 - a11 * a22 * a43 - a12 * a23 * a41 - a13 * a21 * a42);
			b44 = static_cast<d_scalar_t>(a11 * a22 * a33 + a12 * a23 * a31 + a13 * a21 * a32 - a11 * a23 * a32 - a12 * a21 * a33 - a13 * a22 * a31);

			auto detA = determinant44(src_dm16);

			if (maths::is_error(detA) == true) {
				load_identity4(dst_dm16);
				return false;
			}

			detA = static_cast<decltype(detA)>(1.0) / detA;

			for (int i = 0; i < 16; i++) {
				dst_dm16[i] *= static_cast<d_scalar_t>(detA);
			}
			return true;
		}

		////////////////////////////////////////////////
		// 行列の転置
		////////////////////////////////////////////////
		//! @brief IxJ行列を転置しJxI行列に格納する
		//! @param [out] dstJI JxI行列
		//! @param [in] srcIJ IxJ行列
		//! @param [in] I srcIJ行列の行数
		//! @param [in] J srcIJ行列の列数
		//! @return dstJI
		template<typename MatrixD, typename MatrixS>
		inline MatrixD& transpose(MatrixD& dstJI, const MatrixS& srcIJ, const int I, const int J) {
			//i...行
			//j...列
			for (int i = 0; i < I; i++) {
				for (int j = 0; j < J; j++) {
					dstJI[matijI(j, i, J)] = srcIJ[matijI(i, j, I)];
				}
			}
			return dstJI;
		}


		//! @brief 16個のスカラーの一次元配列からなる4*4行列を転置する
		//! @param [in,out] src16 データ real_t[16]
		//! @return src16と同じもの
		template<typename Matrix>
		inline Matrix& transpose44(Matrix& src16) {
			//i...行
			//j...列
			for (int i = 1; i < 4; i++) {
				for (int j = i; j > 0; j--) {
					std::swap(src16[i * 4 + i - j], src16[i * 4 + i - 4 * j]);
				}
			}
			return src16;
		}

		////////////////////////////////////////////////
		// 行列の積
		////////////////////////////////////////////////
		//! @brief IxJ行列 × JxK行列 の積
		//! @param [out] C Cik 計算結果格納先
		//! @param [in] A Aijの行列
		//! @param [in] B Bjkの行列
		//! @param [in] I Aの行列の行数
		//! @param [in] J Aの行列の列数
		//! @param [in] K Bの行列の列数
		//! @return C
		template<typename MatrixC, typename MatrixA, typename MatrixB>
		MatrixC& mult_matrixIJK(MatrixC & C, const MatrixA & A, const MatrixB & B, const int I, const int J, const int K) {

			using scalar_t = typename std::remove_reference<decltype(C[0])>::type;


			for (int i = 0; i < I; i++)
				for (int k = 0; k < K; k++) {
					C[matijI(i, k, I)] = static_cast<scalar_t>(0.0);
					for (int j = 0; j < J; j++) {
						C[matijI(i, k, I)] += static_cast<scalar_t>(A[matijI(i, j, I)] * B[matijI(j, k, J)]);
					}
				}

			return C;
		}


		//! @brief 4x4行列の積
		//! @param [out] mr 結果の格納先。4x4行列
		//! @param [in] ma 4x4行列１
		//! @param [in] mb 4x4行列２
		//! @return mr
		template<typename MatrixD, typename MatrixA, typename MatrixB>
		inline MatrixD& mult_matrix44(MatrixD& mr, const MatrixA& ma, const MatrixB& mb) {
			return mult_matrixIJK(mr, ma, mb, 4, 4, 4);
		}

		//! @brief 4x4行列で四次元ベクトルを変換
		//! @param [out] out4 結果を格納する四次元ベクトル
		//! @param [in] matrix 4x4行列
		//! @param [in] in4 変換する四次元ベクトル
		//! @return out4
		template<typename VectorD4, typename Matrix4, typename VectorS4>
		inline VectorD4& mult_m4_v4(VectorD4& out4, const Matrix4& matrix, const VectorS4& in4)
		{
			int i;

			for (i = 0; i < 4; i++) {
				out4[i] =
					in4[0] * matrix[0 * 4 + i] +
					in4[1] * matrix[1 * 4 + i] +
					in4[2] * matrix[2 * 4 + i] +
					in4[3] * matrix[3 * 4 + i];
			}
			return out4;
		}
	}

}