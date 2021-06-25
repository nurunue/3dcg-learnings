//! @file transformations.hpp
//! @brief OpenGLの関数の実装例
//! @author suzulang.com
//! @date 令和三年六月二十五日

#pragma once

#include "vmmath.hpp"

namespace szl {
	namespace cg3 {

		//! @brief 回転行列作成
		//! @param [out] m 結果を格納する要素数16の配列
		//! @param [in] angle_degree 回転角を度で指定
		//! @param [in] axis_x 回転軸のX成分
		//! @param [in] axis_y 回転軸のY成分
		//! @param [in] axis_z 回転軸のZ成分
		//! @return なし
		template<typename Matrix>
		void rotate_matrix(
			Matrix& m,
			const double angle_degree,
			const double axis_x,
			const double axis_y,
			const double axis_z) {

			using elem_t = typename std::remove_reference<decltype(m[0])>::type;

			double x = axis_x;
			double y = axis_y;
			double z = axis_z;

			//len(x y z) != 0ならnormalizeする
			double len = sqrt(x * x + y * y + z * z);

			if (maths::is_error(len) == false) {
				x = x / len;
				y = y / len;
				z = z / len;
			}
			else {
				mathm::load_identity4(m);
			}


			auto angle_rad = maths::to_radian(angle_degree);

			auto c = std::cos(angle_rad);
			auto s = std::sin(angle_rad);

			m[0] = static_cast<elem_t>(x * x * (1 - c) + c);
			m[1] = static_cast<elem_t>(y * x * (1 - c) + z * s);
			m[2] = static_cast<elem_t>(x * z * (1 - c) - y * s);
			m[3] = static_cast<elem_t>(0.0);

			m[4] = static_cast<elem_t>(x * y * (1 - c) - z * s);
			m[5] = static_cast<elem_t>(y * y * (1 - c) + c);
			m[6] = static_cast<elem_t>(y * z * (1 - c) + x * s);
			m[7] = static_cast<elem_t>(0.0);


			m[8] = static_cast<elem_t>(x * z * (1 - c) + y * s);
			m[9] = static_cast<elem_t>(y * z * (1 - c) - x * s);
			m[10] = static_cast<elem_t>(z * z * (1 - c) + c);
			m[11] = static_cast<elem_t>(0.0);


			m[12] = static_cast<elem_t>(0.0);
			m[13] = static_cast<elem_t>(0.0);
			m[14] = static_cast<elem_t>(0.0);
			m[15] = static_cast<elem_t>(1.0);

		}


		//! @brief 回転行列作成
		//! @param [out] m 結果を格納する要素数16の配列
		//! @param [in] angle_degree 回転角を度で指定
		//! @param [in] axis3 回転軸のX成分
		//! @return なし
		template<typename Matrix, typename Axis, typename ScalarT = double>
		void rotate_matrix(
			Matrix& m,
			const ScalarT angle_degree,
			const Axis& axis3) {

			rotate_matrix(m, angle_degree, axis3[0], axis3[1], axis3[2]);

		}

		//! @brief 移動行列作成
		//! @param [out] m 結果を格納する要素数16の配列
		//! @param [in] x 移動量X成分
		//! @param [in] y 移動量Y成分
		//! @param [in] z 移動量Z成分
		//! @return なし
		template<typename Matrix>
		void translate_matrix(
			Matrix& m,
			const double x,
			const double y,
			const double z
		) {
			using elem_t = typename std::remove_reference<decltype(m[0])>::type;

			m[0] = static_cast<elem_t>(1.0);
			m[1] = static_cast<elem_t>(0.0);
			m[2] = static_cast<elem_t>(0.0);
			m[3] = static_cast<elem_t>(0.0);

			m[4] = static_cast<elem_t>(0.0);
			m[5] = static_cast<elem_t>(1.0);
			m[6] = static_cast<elem_t>(0.0);
			m[7] = static_cast<elem_t>(0.0);

			m[8] = static_cast<elem_t>(0.0);
			m[9] = static_cast<elem_t>(0.0);
			m[10] = static_cast<elem_t>(1.0);
			m[11] = static_cast<elem_t>(0.0);

			m[12] = static_cast<elem_t>(x);
			m[13] = static_cast<elem_t>(y);
			m[14] = static_cast<elem_t>(z);
			m[15] = static_cast<elem_t>(1.0);

		}

		//! @brief 移動行列作成
		//! @param [out] m 結果を格納する要素数16の配列
		//! @param [in] vec3 移動量
		//! @return なし
		template<typename Matrix, typename VectorT>
		void translate_matrix(
			Matrix& m,
			const VectorT& vec3
		) {
			translate_matrix(m, vec3[0], vec3[1], vec3[2]);
		}

		//! @brief 拡大・縮小行列作成
		//! @param [out] m 結果を格納する要素数16の配列
		//! @param [in] sx スケーリング量X
		//! @param [in] sy スケーリング量Y
		//! @param [in] sz スケーリング量Z
		//! @return なし
		template<typename Matrix>
		void scale_matrix(
			Matrix& m,
			const double sx,
			const double sy,
			const double sz
		) {
			using elem_t = typename std::remove_reference<decltype(m[0])>::type;

			m[0] = static_cast<elem_t>(sx);
			m[1] = static_cast<elem_t>(0.0);
			m[2] = static_cast<elem_t>(0.0);
			m[3] = static_cast<elem_t>(0.0);

			m[4] = static_cast<elem_t>(0.0);
			m[5] = static_cast<elem_t>(sy);
			m[6] = static_cast<elem_t>(0.0);
			m[7] = static_cast<elem_t>(0.0);

			m[8] = static_cast<elem_t>(0.0);
			m[9] = static_cast<elem_t>(0.0);
			m[10] = static_cast<elem_t>(sz);
			m[11] = static_cast<elem_t>(0.0);

			m[12] = static_cast<elem_t>(0.0);
			m[13] = static_cast<elem_t>(0.0);
			m[14] = static_cast<elem_t>(0.0);
			m[15] = static_cast<elem_t>(1.0);

		}
		//! @brief 拡大・縮小行列作成
		//! @param [out] m 結果を格納する要素数16の配列
		//! @param [in] vec3 スケーリング量
		//! @return なし
		template<typename Matrix, typename VectorT>
		void scale_matrix(
			Matrix& m,
			const VectorT& vec3
		) {
			scale_matrix(m, vec3[0], vec3[1], vec3[2]);
		}





		//////////////////////////////////////////////////

		//! @brief 透視投影行列の作成 (視野角で指定)
		//! @param [out] m 結果を格納する要素数16の配列
		//! @param [in] fovy_degree 視野角
		//! @param [in] aspect アスペクト比
		//! @param [in] zNear 一番近いz位置
		//! @param [in] zFar 一番遠いz位置
		//! @return なし
		template<typename Matrix>
		void perspective_matrix(
			Matrix& m,
			double fovy_degree,
			double aspect,
			double zNear,
			double zFar) {

			using elem_t = typename std::remove_reference<decltype(m[0])>::type;


			double fovy_rad = maths::to_radian(fovy_degree);


			double f = maths::cot(fovy_rad / 2.0);

			m[0] = static_cast<elem_t>(f / aspect);
			m[1] = static_cast<elem_t>(0.0);
			m[2] = static_cast<elem_t>(0.0);
			m[3] = static_cast<elem_t>(0.0);

			m[4] = static_cast<elem_t>(0.0);
			m[5] = static_cast<elem_t>(f);
			m[6] = static_cast<elem_t>(0.0);
			m[7] = static_cast<elem_t>(0.0);

			m[8] = static_cast<elem_t>(0.0);
			m[9] = static_cast<elem_t>(0.0);
			m[10] = static_cast<elem_t>((zFar + zNear) / (zNear - zFar));
			m[11] = static_cast<elem_t>(-1.0);

			m[12] = static_cast<elem_t>(0.0);
			m[13] = static_cast<elem_t>(0.0);
			m[14] = static_cast<elem_t>((2 * zFar * zNear) / (zNear - zFar));
			m[15] = static_cast<elem_t>(0.0);

		}

		//! @brief 平行投影行列作成
		//! @param [out] m 結果を格納する要素数16の配列
		//! @param [in] left 左の垂直方向のクリッピング平面の座標
		//! @param [in] right 垂直方向のクリッピング平面の座標
		//! @param [in] bottom 下方向のクリッピング平面の座標
		//! @param [in] top 上部の水平方向のクリッピング計画の座標
		//! @param [in] znear 最も近い深度のクリッピング平面への距離
		//! @param [in] zfar より深い深度のクリッピング平面への距離
		//! @return
		template<typename Matrix>
		void ortho_matrix(
			Matrix& m,
			const double left,
			const double right,
			const double bottom,
			const double top,
			const double znear,
			const double zfar) {

			using elem_t = typename std::remove_reference<decltype(m[0])>::type;

			double tx = -(right + left) / (right - left);
			double ty = -(top + bottom) / (top - bottom);
			double tz = -(zfar + znear) / (zfar - znear);

			//infが出たときは単位行列を返す模様
			if (std::isinf(tx) || std::isinf(ty) || std::isinf(tz)) {
				mathm::load_identity4(m);

				return;
			}



			m[0] = static_cast<elem_t>(2.0 / (right - left));
			m[1] = static_cast<elem_t>(0.0);
			m[2] = static_cast<elem_t>(0.0);
			m[3] = static_cast<elem_t>(0.0);

			m[4] = static_cast<elem_t>(0.0);
			m[5] = static_cast<elem_t>(2.0 / (top - bottom));
			m[6] = static_cast<elem_t>(0.0);
			m[7] = static_cast<elem_t>(0.0);

			m[8] = static_cast<elem_t>(0.0);
			m[9] = static_cast<elem_t>(0.0);
			m[10] = static_cast<elem_t>(-2.0 / (zfar - znear));
			m[11] = static_cast<elem_t>(0.0);

			m[12] = static_cast<elem_t>(tx);
			m[13] = static_cast<elem_t>(ty);
			m[14] = static_cast<elem_t>(tz);
			m[15] = static_cast<elem_t>(1.0);

		}


		//! @brief gluLookAtの自前実装版
		//! @param [out] M 結果の格納先
		//! @param [in] eyeX カメラ位置X
		//! @param [in] eyeY カメラ位置Y
		//! @param [in] eyeZ カメラ位置Z
		//! @param [in] centerX 注視点X
		//! @param [in] centerY 注視点Y
		//! @param [in] centerZ 注視点Z
		//! @param [in] upX 上方向X
		//! @param [in] upY 上方向Y
		//! @param [in] upZ 上方向Z
		//! @return なし
		template<typename Matrix>
		void lookat_matrix(
			Matrix& M,
			double eyeX,
			double eyeY,
			double eyeZ,
			double centerX,
			double centerY,
			double centerZ,
			double upX,
			double upY,
			double upZ)
		{
			double f[3] = {
			  centerX - eyeX,
			  centerY - eyeY,
			  centerZ - eyeZ
			};
			mathv::normalize3(f);

			double UP[3] = { upX,upY,upZ };
			mathv::normalize3(UP);

			double s[3];
			mathv::outer3(s, f, UP);
			mathv::normalize3(s);

			double u[3];
			mathv::outer3(u, s, f);

			double M0[16];
			M0[0] = s[0];
			M0[1] = u[0];
			M0[2] = -f[0];
			M0[3] = 0;

			M0[4] = s[1];
			M0[5] = u[1];
			M0[6] = -f[1];
			M0[7] = 0;

			M0[8] = s[2];
			M0[9] = u[2];
			M0[10] = -f[2];
			M0[11] = 0;

			M0[12] = 0;
			M0[13] = 0;
			M0[14] = 0;
			M0[15] = 1;


			/////////////////////////////
			/////////////////////////////
			double E[16];

			translate_matrix(E, -eyeX, -eyeY, -eyeZ);

			mathm::mult_matrix44(M, M0, E);

			if (mathv::chek_array_nan(M, 16) == true) {
				mathm::load_identity4(M);
			}

		}

		//! @brief gluLookAtの自前実装版
		//! @param [out] M 結果の格納先
		//! @param [in] eye3 カメラ位置
		//! @param [in] center3 注視点
		//! @param [in] up3 上方向
		//! @return なし
		template<typename Matrix, typename VectorE, typename VectorC, typename VectorU>
		void lookat_matrix(
			Matrix& M,
			const VectorE& eye3,
			const VectorC& center3,
			const VectorU& up3
		)
		{
			lookat_matrix(M,
				eye3[0], eye3[1], eye3[2],
				center3[0], center3[1], center3[2],
				up3[0], up3[1], up3[2]
			);
		}

	}
}