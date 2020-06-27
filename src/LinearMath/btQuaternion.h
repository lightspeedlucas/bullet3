/*
Copyright (c) 2003-2006 Gino van den Bergen / Erwin Coumans  http://continuousphysics.com/Bullet/

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose, 
including commercial applications, and to alter it and redistribute it freely, 
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/

#ifndef BT_SIMD__QUATERNION_H_
#define BT_SIMD__QUATERNION_H_

#include "btVector3.h"
#include "btQuadWord.h"

#ifdef BT_USE_DOUBLE_PRECISION
#define btQuaternionData btQuaternionDoubleData
#define btQuaternionDataName "btQuaternionDoubleData"
#else
#define btQuaternionData btQuaternionFloatData
#define btQuaternionDataName "btQuaternionFloatData"
#endif  //BT_USE_DOUBLE_PRECISION

#ifdef BT_USE_SSE

//const __m128 ATTRIBUTE_ALIGNED16(vOnes) = {1.0f, 1.0f, 1.0f, 1.0f};
#define vOnes (_mm_set_ps(1.0f, 1.0f, 1.0f, 1.0f))

#endif

#if defined(BT_USE_SSE)

#define vQInv (_mm_set_ps(+0.0f, -0.0f, -0.0f, -0.0f))
#define vPPPM (_mm_set_ps(-0.0f, +0.0f, +0.0f, +0.0f))

#elif defined(BT_USE_NEON)

const btSimdFloat4 ATTRIBUTE_ALIGNED16(vQInv) = {-0.0f, -0.0f, -0.0f, +0.0f};
const btSimdFloat4 ATTRIBUTE_ALIGNED16(vPPPM) = {+0.0f, +0.0f, +0.0f, -0.0f};

#endif

/**@brief The btQuaternion implements quaternion to perform linear algebra rotations in combination with btMatrix3x3, btVector3 and btTransform. */
class btQuaternion : public btQuadWord
{
public:
	/**@brief No initialization constructor */
	btQuaternion() {}

	//		template <typename btScalar>
	//		explicit Quaternion(const btScalar *v) : Tuple4<btScalar>(v) {}
	/**@brief Constructor from scalars */
	btQuaternion(const btScalar& _x, const btScalar& _y, const btScalar& _z, const btScalar& _w)
		: btQuadWord(_x, _y, _z, _w)
	{
	}
	/**@brief Axis angle Constructor
   * @param axis The axis which the rotation is around
   * @param angle The magnitude of the rotation around the angle (Radians) */
	btQuaternion(const btVector3& _axis, const btScalar& _angle)
	{
		setRotation(_axis, _angle);
	}
	/**@brief Constructor from Euler angles
   * @param yaw Angle around Y unless BT_EULER_DEFAULT_ZYX defined then Z
   * @param pitch Angle around X unless BT_EULER_DEFAULT_ZYX defined then Y
   * @param roll Angle around Z unless BT_EULER_DEFAULT_ZYX defined then X */
	btQuaternion(const btScalar& yaw, const btScalar& pitch, const btScalar& roll)
	{
#ifndef BT_EULER_DEFAULT_ZYX
		setEuler(yaw, pitch, roll);
#else
		setEulerZYX(yaw, pitch, roll);
#endif
	}

	// Needs to be declared explicitly now (as of fixed point implementation)
	// for some reason
	btQuaternion(const btQuaternion& q) : btQuadWord(q.x(), q.y(), q.z(), q.w())
	{
	}
	/**@brief Set the rotation using axis angle notation 
   * @param axis The axis around which to rotate
   * @param angle The magnitude of the rotation in Radians */
	void setRotation(const btVector3& axis, const btScalar& _angle)
	{
		btScalar d = axis.length();
		btAssert(d != btScalar(0.0));
		btScalar s = btSin(_angle * btScalar(0.5)) / d;
		setValue(axis.x() * s, axis.y() * s, axis.z() * s,
				 btCos(_angle * btScalar(0.5)));
	}
	/**@brief Set the quaternion using Euler angles
   * @param yaw Angle around Y
   * @param pitch Angle around X
   * @param roll Angle around Z */
	void setEuler(const btScalar& yaw, const btScalar& pitch, const btScalar& roll)
	{
		btScalar halfYaw = btScalar(yaw) * btScalar(0.5);
		btScalar halfPitch = btScalar(pitch) * btScalar(0.5);
		btScalar halfRoll = btScalar(roll) * btScalar(0.5);
		btScalar cosYaw = btCos(halfYaw);
		btScalar sinYaw = btSin(halfYaw);
		btScalar cosPitch = btCos(halfPitch);
		btScalar sinPitch = btSin(halfPitch);
		btScalar cosRoll = btCos(halfRoll);
		btScalar sinRoll = btSin(halfRoll);
		setValue(cosRoll * sinPitch * cosYaw + sinRoll * cosPitch * sinYaw,
				 cosRoll * cosPitch * sinYaw - sinRoll * sinPitch * cosYaw,
				 sinRoll * cosPitch * cosYaw - cosRoll * sinPitch * sinYaw,
				 cosRoll * cosPitch * cosYaw + sinRoll * sinPitch * sinYaw);
	}
	/**@brief Set the quaternion using euler angles 
   * @param yaw Angle around Z
   * @param pitch Angle around Y
   * @param roll Angle around X */
	void setEulerZYX(const btScalar& yawZ, const btScalar& pitchY, const btScalar& rollX)
	{
		btScalar halfYaw = btScalar(yawZ) * btScalar(0.5);
		btScalar halfPitch = btScalar(pitchY) * btScalar(0.5);
		btScalar halfRoll = btScalar(rollX) * btScalar(0.5);
		btScalar cosYaw = btCos(halfYaw);
		btScalar sinYaw = btSin(halfYaw);
		btScalar cosPitch = btCos(halfPitch);
		btScalar sinPitch = btSin(halfPitch);
		btScalar cosRoll = btCos(halfRoll);
		btScalar sinRoll = btSin(halfRoll);
		setValue(sinRoll * cosPitch * cosYaw - cosRoll * sinPitch * sinYaw,   //x
				 cosRoll * sinPitch * cosYaw + sinRoll * cosPitch * sinYaw,   //y
				 cosRoll * cosPitch * sinYaw - sinRoll * sinPitch * cosYaw,   //z
				 cosRoll * cosPitch * cosYaw + sinRoll * sinPitch * sinYaw);  //formerly yzx
	}

	/**@brief Get the euler angles from this quaternion
	   * @param yaw Angle around Z
	   * @param pitch Angle around Y
	   * @param roll Angle around X */
	void getEulerZYX(btScalar& yawZ, btScalar& pitchY, btScalar& rollX) const
	{
		btScalar squ;
		btScalar sqx;
		btScalar sqy;
		btScalar sqz;
		btScalar sarg;
		sqx = m_floats[0] * m_floats[0];
		sqy = m_floats[1] * m_floats[1];
		sqz = m_floats[2] * m_floats[2];
		squ = m_floats[3] * m_floats[3];
		sarg = btScalar(-2.) * (m_floats[0] * m_floats[2] - m_floats[3] * m_floats[1]);

		// If the pitch angle is PI/2 or -PI/2, we can only compute
		// the sum roll + yaw.  However, any combination that gives
		// the right sum will produce the correct orientation, so we
		// set rollX = 0 and compute yawZ.
		if (sarg <= -btScalar(0.99999))
		{
			pitchY = btScalar(-0.5) * SIMD_PI;
			rollX = 0;
			yawZ = btScalar(2) * btAtan2(m_floats[0], -m_floats[1]);
		}
		else if (sarg >= btScalar(0.99999))
		{
			pitchY = btScalar(0.5) * SIMD_PI;
			rollX = 0;
			yawZ = btScalar(2) * btAtan2(-m_floats[0], m_floats[1]);
		}
		else
		{
			pitchY = btAsin(sarg);
			rollX = btAtan2(btScalar(2) * (m_floats[1] * m_floats[2] + m_floats[3] * m_floats[0]), squ - sqx - sqy + sqz);
			yawZ = btAtan2(btScalar(2) * (m_floats[0] * m_floats[1] + m_floats[3] * m_floats[2]), squ + sqx - sqy - sqz);
		}
	}

	/**@brief Add two quaternions
   * @param q The quaternion to add to this one */
	SIMD_FORCE_INLINE btQuaternion& operator+=(const btQuaternion& q)
	{
		m_floats[0] += q.x();
		m_floats[1] += q.y();
		m_floats[2] += q.z();
		m_floats[3] += q.m_floats[3];
		return *this;
	}

	/**@brief Subtract out a quaternion
   * @param q The quaternion to subtract from this one */
	btQuaternion& operator-=(const btQuaternion& q)
	{
		m_floats[0] -= q.x();
		m_floats[1] -= q.y();
		m_floats[2] -= q.z();
		m_floats[3] -= q.m_floats[3];
		return *this;
	}

	/**@brief Scale this quaternion
   * @param s The scalar to scale by */
	btQuaternion& operator*=(const btScalar& s)
	{
		m_floats[0] *= s;
		m_floats[1] *= s;
		m_floats[2] *= s;
		m_floats[3] *= s;
		return *this;
	}

	/**@brief Multiply this quaternion by q on the right
   * @param q The other quaternion 
   * Equivilant to this = this * q */
	btQuaternion& operator*=(const btQuaternion& q)
	{
		setValue(
			m_floats[3] * q.x() + m_floats[0] * q.m_floats[3] + m_floats[1] * q.z() - m_floats[2] * q.y(),
			m_floats[3] * q.y() + m_floats[1] * q.m_floats[3] + m_floats[2] * q.x() - m_floats[0] * q.z(),
			m_floats[3] * q.z() + m_floats[2] * q.m_floats[3] + m_floats[0] * q.y() - m_floats[1] * q.x(),
			m_floats[3] * q.m_floats[3] - m_floats[0] * q.x() - m_floats[1] * q.y() - m_floats[2] * q.z());
		return *this;
	}
	/**@brief Return the dot product between this quaternion and another
   * @param q The other quaternion */
	btScalar dot(const btQuaternion& q) const
	{
#if defined BT_USE_SIMD_VECTOR3 && defined(BT_USE_SSE_IN_API) && defined(BT_USE_SSE)
		__m128 vd;

		vd = _mm_mul_ps(mVec128, q.mVec128);

		__m128 t = _mm_movehl_ps(vd, vd);
		vd = _mm_add_ps(vd, t);
		t = _mm_shuffle_ps(vd, vd, 0x55);
		vd = _mm_add_ss(vd, t);

		return _mm_cvtss_f32(vd);
#elif defined(BT_USE_NEON)
		float32x4_t vd = vmulq_f32(mVec128, q.mVec128);
		float32x2_t x = vpadd_f32(vget_low_f32(vd), vget_high_f32(vd));
		x = vpadd_f32(x, x);
		return vget_lane_f32(x, 0);
#else
		return m_floats[0] * q.x() +
			   m_floats[1] * q.y() +
			   m_floats[2] * q.z() +
			   m_floats[3] * q.m_floats[3];
#endif
	}

	/**@brief Return the length squared of the quaternion */
	btScalar length2() const
	{
		return dot(*this);
	}

	/**@brief Return the length of the quaternion */
	btScalar length() const
	{
		return btSqrt(length2());
	}
	btQuaternion& safeNormalize()
	{
		btScalar l2 = length2();
		if (l2 > btScalar(SIMD_EPSILON))
		{
			normalize();
		}
		return *this;
	}
	/**@brief Normalize the quaternion 
   * Such that x^2 + y^2 + z^2 +w^2 = 1 */
	btQuaternion& normalize()
	{
		return *this /= length();
	}

	/**@brief Return a scaled version of this quaternion
   * @param s The scale factor */
	SIMD_FORCE_INLINE btQuaternion
	operator*(const btScalar& s) const
	{
		return btQuaternion(x() * s, y() * s, z() * s, m_floats[3] * s);
	}

	/**@brief Return an inversely scaled versionof this quaternion
   * @param s The inverse scale factor */
	btQuaternion operator/(const btScalar& s) const
	{
		btAssert(s != btScalar(0.0));
		return *this * (btScalar(1.0) / s);
	}

	/**@brief Inversely scale this quaternion
   * @param s The scale factor */
	btQuaternion& operator/=(const btScalar& s)
	{
		btAssert(s != btScalar(0.0));
		return *this *= btScalar(1.0) / s;
	}

	/**@brief Return a normalized version of this quaternion */
	btQuaternion normalized() const
	{
		return *this / length();
	}
	/**@brief Return the ***half*** angle between this quaternion and the other
   * @param q The other quaternion */
	btScalar angle(const btQuaternion& q) const
	{
		btScalar s = btSqrt(length2() * q.length2());
		btAssert(s != btScalar(0.0));
		return btAcos(dot(q) / s);
	}

	/**@brief Return the angle between this quaternion and the other along the shortest path
	* @param q The other quaternion */
	btScalar angleShortestPath(const btQuaternion& q) const
	{
		btScalar s = btSqrt(length2() * q.length2());
		btAssert(s != btScalar(0.0));
		if (dot(q) < btScalar(0))  // Take care of long angle case see http://en.wikipedia.org/wiki/Slerp
			return btAcos(dot(-q) / s) * btScalar(2.0);
		else
			return btAcos(dot(q) / s) * btScalar(2.0);
	}

	/**@brief Return the angle [0, 2Pi] of rotation represented by this quaternion */
	btScalar getAngle() const
	{
		btScalar s = btScalar(2.) * btAcos(m_floats[3]);
		return s;
	}

	/**@brief Return the angle [0, Pi] of rotation represented by this quaternion along the shortest path */
	btScalar getAngleShortestPath() const
	{
		btScalar s;
		if (m_floats[3] >= btScalar(0))
			s = btScalar(2.) * btAcos(m_floats[3]);
		else
			s = btScalar(2.) * btAcos(-m_floats[3]);
		return s;
	}

	/**@brief Return the axis of the rotation represented by this quaternion */
	btVector3 getAxis() const
	{
		btScalar s_squared = btScalar(1) - m_floats[3] * m_floats[3];

		if (s_squared < btScalar(10.) * btScalar(SIMD_EPSILON))  //Check for divide by zero
			return btVector3(1.0, 0.0, 0.0);           // Arbitrary
		btScalar s = btScalar(1) / btSqrt(s_squared);
		return btVector3(m_floats[0] * s, m_floats[1] * s, m_floats[2] * s);
	}

	/**@brief Return the inverse of this quaternion */
	btQuaternion inverse() const
	{
		return btQuaternion(-m_floats[0], -m_floats[1], -m_floats[2], m_floats[3]);
	}

	/**@brief Return the sum of this quaternion and the other 
   * @param q2 The other quaternion */
	SIMD_FORCE_INLINE btQuaternion
	operator+(const btQuaternion& q2) const
	{
		const btQuaternion& q1 = *this;
		return btQuaternion(q1.x() + q2.x(), q1.y() + q2.y(), q1.z() + q2.z(), q1.m_floats[3] + q2.m_floats[3]);
	}

	/**@brief Return the difference between this quaternion and the other 
   * @param q2 The other quaternion */
	SIMD_FORCE_INLINE btQuaternion
	operator-(const btQuaternion& q2) const
	{
		const btQuaternion& q1 = *this;
		return btQuaternion(q1.x() - q2.x(), q1.y() - q2.y(), q1.z() - q2.z(), q1.m_floats[3] - q2.m_floats[3]);
	}

	/**@brief Return the negative of this quaternion 
   * This simply negates each element */
	SIMD_FORCE_INLINE btQuaternion operator-() const
	{
		const btQuaternion& q2 = *this;
		return btQuaternion(-q2.x(), -q2.y(), -q2.z(), -q2.m_floats[3]);
	}
	/**@todo document this and it's use */
	SIMD_FORCE_INLINE btQuaternion farthest(const btQuaternion& qd) const
	{
		btQuaternion diff, sum;
		diff = *this - qd;
		sum = *this + qd;
		if (diff.dot(diff) > sum.dot(sum))
			return qd;
		return (-qd);
	}

	/**@todo document this and it's use */
	SIMD_FORCE_INLINE btQuaternion nearest(const btQuaternion& qd) const
	{
		btQuaternion diff, sum;
		diff = *this - qd;
		sum = *this + qd;
		if (diff.dot(diff) < sum.dot(sum))
			return qd;
		return (-qd);
	}

	/**@brief Return the quaternion which is the result of Spherical Linear Interpolation between this and the other quaternion
   * @param q The other quaternion to interpolate with 
   * @param t The ratio between this and q to interpolate.  If t = 0 the result is this, if t=1 the result is q.
   * Slerp interpolates assuming constant velocity.  */
	btQuaternion slerp(const btQuaternion& q, const btScalar& t) const
	{
		const btScalar magnitude = btSqrt(length2() * q.length2());
		btAssert(magnitude > btScalar(0));

		const btScalar product = dot(q) / magnitude;
		const btScalar absproduct = btFabs(product);

		if (absproduct < btScalar(1.0 - SIMD_EPSILON))
		{
			// Take care of long angle case see http://en.wikipedia.org/wiki/Slerp
			const btScalar theta = btAcos(absproduct);
			const btScalar d = btSin(theta);
			btAssert(d > btScalar(0));

			const btScalar sign = (product < btScalar(0)) ? btScalar(-1) : btScalar(1);
			const btScalar s0 = btSin((btScalar(1.0) - t) * theta) / d;
			const btScalar s1 = btSin(sign * t * theta) / d;

			return btQuaternion(
				(m_floats[0] * s0 + q.x() * s1),
				(m_floats[1] * s0 + q.y() * s1),
				(m_floats[2] * s0 + q.z() * s1),
				(m_floats[3] * s0 + q.w() * s1));
		}
		else
		{
			return *this;
		}
	}

	static const btQuaternion& getIdentity()
	{
		static const btQuaternion identityQuat(btScalar(0.), btScalar(0.), btScalar(0.), btScalar(1.));
		return identityQuat;
	}

	SIMD_FORCE_INLINE const btScalar& getW() const { return m_floats[3]; }

	SIMD_FORCE_INLINE void serialize(struct btQuaternionData& dataOut) const;

	SIMD_FORCE_INLINE void deSerialize(const struct btQuaternionFloatData& dataIn);

	SIMD_FORCE_INLINE void deSerialize(const struct btQuaternionDoubleData& dataIn);

	SIMD_FORCE_INLINE void serializeFloat(struct btQuaternionFloatData& dataOut) const;

	SIMD_FORCE_INLINE void deSerializeFloat(const struct btQuaternionFloatData& dataIn);

	SIMD_FORCE_INLINE void serializeDouble(struct btQuaternionDoubleData& dataOut) const;

	SIMD_FORCE_INLINE void deSerializeDouble(const struct btQuaternionDoubleData& dataIn);
};

/**@brief Return the product of two quaternions */
SIMD_FORCE_INLINE btQuaternion
operator*(const btQuaternion& q1, const btQuaternion& q2)
{
	return btQuaternion(
		q1.w() * q2.x() + q1.x() * q2.w() + q1.y() * q2.z() - q1.z() * q2.y(),
		q1.w() * q2.y() + q1.y() * q2.w() + q1.z() * q2.x() - q1.x() * q2.z(),
		q1.w() * q2.z() + q1.z() * q2.w() + q1.x() * q2.y() - q1.y() * q2.x(),
		q1.w() * q2.w() - q1.x() * q2.x() - q1.y() * q2.y() - q1.z() * q2.z());
}

SIMD_FORCE_INLINE btQuaternion
operator*(const btQuaternion& q, const btVector3& w)
{
	return btQuaternion(
		q.w() * w.x() + q.y() * w.z() - q.z() * w.y(),
		q.w() * w.y() + q.z() * w.x() - q.x() * w.z(),
		q.w() * w.z() + q.x() * w.y() - q.y() * w.x(),
		-q.x() * w.x() - q.y() * w.y() - q.z() * w.z());
}

SIMD_FORCE_INLINE btQuaternion
operator*(const btVector3& w, const btQuaternion& q)
{
	return btQuaternion(
		+w.x() * q.w() + w.y() * q.z() - w.z() * q.y(),
		+w.y() * q.w() + w.z() * q.x() - w.x() * q.z(),
		+w.z() * q.w() + w.x() * q.y() - w.y() * q.x(),
		-w.x() * q.x() - w.y() * q.y() - w.z() * q.z());
}

/**@brief Calculate the dot product between two quaternions */
SIMD_FORCE_INLINE btScalar
dot(const btQuaternion& q1, const btQuaternion& q2)
{
	return q1.dot(q2);
}

/**@brief Return the length of a quaternion */
SIMD_FORCE_INLINE btScalar
length(const btQuaternion& q)
{
	return q.length();
}

/**@brief Return the angle between two quaternions*/
SIMD_FORCE_INLINE btScalar
btAngle(const btQuaternion& q1, const btQuaternion& q2)
{
	return q1.angle(q2);
}

/**@brief Return the inverse of a quaternion*/
SIMD_FORCE_INLINE btQuaternion
inverse(const btQuaternion& q)
{
	return q.inverse();
}

/**@brief Return the result of spherical linear interpolation betwen two quaternions 
 * @param q1 The first quaternion
 * @param q2 The second quaternion 
 * @param t The ration between q1 and q2.  t = 0 return q1, t=1 returns q2 
 * Slerp assumes constant velocity between positions. */
SIMD_FORCE_INLINE btQuaternion
slerp(const btQuaternion& q1, const btQuaternion& q2, const btScalar& t)
{
	return q1.slerp(q2, t);
}

SIMD_FORCE_INLINE btVector3
quatRotate(const btQuaternion& rotation, const btVector3& v)
{
	btQuaternion q = rotation * v;
	q *= rotation.inverse();
	return btVector3(q.getX(), q.getY(), q.getZ());
}

SIMD_FORCE_INLINE btQuaternion
shortestArcQuat(const btVector3& v0, const btVector3& v1)  // Game Programming Gems 2.10. make sure v0,v1 are normalized
{
	btVector3 c = v0.cross(v1);
	btScalar d = v0.dot(v1);

	if (d < btScalar(-1) + btScalar(SIMD_EPSILON))
	{
		btVector3 n, unused;
		btPlaneSpace1(v0, n, unused);
		return btQuaternion(n.x(), n.y(), n.z(), 0.0f);  // just pick any vector that is orthogonal to v0
	}

	btScalar s = btSqrt((btScalar(1) + d) * btScalar(2));
	btScalar rs = btScalar(1) / s;

	return btQuaternion(c.getX() * rs, c.getY() * rs, c.getZ() * rs, s * btScalar(0.5));
}

SIMD_FORCE_INLINE btQuaternion
shortestArcQuatNormalize2(btVector3& v0, btVector3& v1)
{
	v0.normalize();
	v1.normalize();
	return shortestArcQuat(v0, v1);
}

struct btQuaternionFloatData
{
	float m_floats[4];
};

struct btQuaternionDoubleData
{
	double m_floats[4];
};

SIMD_FORCE_INLINE void btQuaternion::serializeFloat(struct btQuaternionFloatData& dataOut) const
{
	///could also do a memcpy, check if it is worth it
	for (int i = 0; i < 4; i++)
		dataOut.m_floats[i] = float(m_floats[i].ToFloat());
}

SIMD_FORCE_INLINE void btQuaternion::deSerializeFloat(const struct btQuaternionFloatData& dataIn)
{
	for (int i = 0; i < 4; i++)
		m_floats[i] = btScalar(dataIn.m_floats[i]);
}

SIMD_FORCE_INLINE void btQuaternion::serializeDouble(struct btQuaternionDoubleData& dataOut) const
{
	///could also do a memcpy, check if it is worth it
	for (int i = 0; i < 4; i++)
		dataOut.m_floats[i] = double(m_floats[i].ToFloat());
}

SIMD_FORCE_INLINE void btQuaternion::deSerializeDouble(const struct btQuaternionDoubleData& dataIn)
{
	for (int i = 0; i < 4; i++)
		m_floats[i] = btScalar(dataIn.m_floats[i]);
}

SIMD_FORCE_INLINE void btQuaternion::serialize(struct btQuaternionData& dataOut) const
{
	///could also do a memcpy, check if it is worth it
	for (int i = 0; i < 4; i++)
		dataOut.m_floats[i] = m_floats[i].ToFloat();
}

SIMD_FORCE_INLINE void btQuaternion::deSerialize(const struct btQuaternionFloatData& dataIn)
{
	for (int i = 0; i < 4; i++)
		m_floats[i] = (btScalar)dataIn.m_floats[i];
}

SIMD_FORCE_INLINE void btQuaternion::deSerialize(const struct btQuaternionDoubleData& dataIn)
{
	for (int i = 0; i < 4; i++)
		m_floats[i] = (btScalar)dataIn.m_floats[i];
}

#endif  //BT_SIMD__QUATERNION_H_
