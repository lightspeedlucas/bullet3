#pragma once

#include <stdint.h>
#define _USE_MATH_DEFINES
#include <math.h>

class Q
{
	typedef int64_t BaseType;
	typedef uint64_t UnsignedBaseType;
	typedef uint32_t UnsignedHalfType;

	BaseType Value;

public:
	static constexpr size_t FractionalBits = 16;
	static constexpr size_t Scale = (1 << FractionalBits);
	static constexpr UnsignedBaseType FractionalMask = (1ULL << FractionalBits) - 1;
	static constexpr UnsignedBaseType IntegerMask = ~FractionalMask;
	static constexpr float Precision = 1.0f / (float) Scale;

	// Constructors
	Q() : Value(0)
	{
	}

	Q(int16_t I)
	{
		BaseType Temp = I;
		Value = Temp << FractionalBits;
	}

	Q(uint16_t I)
	{
		BaseType Temp = I;
		Value = Temp << FractionalBits;
	}

	Q(int32_t I)
	{
		BaseType Temp = I;
		Value = Temp << FractionalBits;
	}

	Q(uint32_t I)
	{
		BaseType Temp = I;
		Value = Temp << FractionalBits;
	}

	Q(BaseType I)
	{
		Value = I << FractionalBits;
	}

	Q(UnsignedBaseType I)
	{
		Value = I << FractionalBits;
	}

	Q(float F)
	{
		Value = (BaseType) (F * Scale);
	}

	Q(double D)
	{
		Value = (BaseType) (D * Scale);
	}

	Q(const Q& S)
	{
		Value = S.Value;
	}

	static Q FromRaw(BaseType RawValue)
	{
		Q Res;
		Res.Value = RawValue;
		return Res;
	}

	// Typecasts
	operator short() const
	{
		return (short)ToInt();
	}

	operator unsigned short() const
	{
		return (unsigned short)ToInt();
	}

	operator int() const
	{
		return (int)ToInt();
	}

	operator unsigned int() const
	{
		return (unsigned int)ToInt();
	}

	operator float() const
	{
		return ToFloat();
	}

	// Arithmetic operations
	Q operator+(Q Other) const
	{
		Q Res;
		Res.Value = Value + Other.Value;
		return Res;
	}

	Q operator-(Q Other) const
	{
		Q Res;
		Res.Value = Value - Other.Value;
		return Res;
	}

	Q operator-() const
	{
		return Negate(*this);
	}

	Q operator+() const
	{
		return *this;
	}

	// Taken from https://www.drdobbs.com/cpp/optimizing-math-intensive-applications-w/207000448?pgno=1
	Q operator*(Q Other) const
	{
		bool const OtherNegative = Other.Value < 0;
		bool const ThisNegative = Value < 0;
		bool const Negate = OtherNegative ^ ThisNegative;
		UnsignedBaseType const OtherVal = OtherNegative ? -Other.Value : Other.Value;
		UnsignedBaseType const Self = ThisNegative ? -Value : Value;

		BaseType Res;
		if (UnsignedBaseType const SelfUpper = (Self >> 32))
		{
			Res = (SelfUpper * OtherVal) << (32 - FractionalBits);
		}
		else
		{
			Res = 0;
		}

		if (UnsignedBaseType const SelfLower = (Self & 0xffffffff))
		{
			UnsignedHalfType const OtherUpper = (UnsignedHalfType)(OtherVal >> 32);
			UnsignedHalfType const OtherLower = (UnsignedHalfType)(OtherVal & 0xffffffff);
			UnsignedBaseType const LowerSelfUpperOtherRes = SelfLower * OtherUpper;
			UnsignedBaseType const LowerSelfLowerOtherRes = SelfLower * OtherLower;
			Res += (LowerSelfUpperOtherRes << (32 - FractionalBits))
				+ (LowerSelfLowerOtherRes >> FractionalBits);
		}

		if (Negate)
		{
			Res = -Res;
		}
		return FromRaw(Res);
	}

	Q operator/(Q Other) const
	{
		// Division by 0 returns the absolute lowest value possible
		if (Other == Q(0))
		{
			return FromRaw(0x8000000000000000LL);
		}

		// Taken from https://www.drdobbs.com/cpp/optimizing-math-intensive-applications-w/207000448?pgno=1
		bool const NegateThis = (Value < 0);
		bool const NegateDivisor = (Other.Value < 0);
		bool const Negate = NegateThis ^ NegateDivisor;
		UnsignedBaseType a = NegateThis ? -Value : Value;
		UnsignedBaseType b = NegateDivisor ? -Other.Value : Other.Value;

		UnsignedBaseType res = 0;

		UnsignedBaseType temp = b;
		bool const AIsLarge = a > b;
		size_t shift = FractionalBits;

		if (AIsLarge)
		{
			UnsignedBaseType const half_a = a >> 1;
			while (temp < half_a)
			{
				temp <<= 1;
				++shift;
			}
		}
		UnsignedBaseType d = 1I64 << shift;
		if (AIsLarge)
		{
			a -= temp;
			res += d;
		}

		while (a && temp && shift)
		{
			unsigned right_shift = 0;
			while (right_shift < shift && (temp > a))
			{
				temp >>= 1;
				++right_shift;
			}
			d >>= right_shift;
			shift -= right_shift;
			a -= temp;
			res += d;
		}
		return FromRaw(Negate ? -(BaseType)res : res);
	}

	// Assignment operators
	Q& operator+=(const Q& Rhs)
	{
		*this = *this + Rhs;
		return *this;
	}

	Q& operator-=(const Q& Rhs)
	{
		*this = *this - Rhs;
		return *this;
	}

	Q& operator*=(const Q& Rhs)
	{
		*this = *this * Rhs;
		return *this;
	}

	Q& operator/=(const Q& Rhs)
	{
		*this = *this / Rhs;
		return *this;
	}

	// Comparison operators
	bool operator==(Q Other) const
	{
		return Value == Other.Value;
	}

	bool operator!=(Q Other) const
	{
		return Value != Other.Value;
	}

	bool operator>(Q Other) const
	{
		return Value > Other.Value;
	}

	bool operator>=(Q Other) const
	{
		return Value >= Other.Value;
	}

	bool operator<(Q Other) const
	{
		return Value < Other.Value;
	}

	bool operator<=(Q Other) const
	{
		return Value <= Other.Value;
	}

	// Convenient comparisons
	bool operator==(int Other) const
	{
		return *this == Q(Other);
	}

	bool operator==(float Other) const
	{
		return *this == Q(Other);
	}

	bool operator!=(int Other) const
	{
		return *this != Q(Other);
	}

	bool operator!=(float Other) const
	{
		return *this != Q(Other);
	}

	bool operator>(int Other) const
	{
		return *this > Q(Other);
	}

	bool operator>(float Other) const
	{
		return *this > Q(Other);
	}

	bool operator>=(int Other) const
	{
		return *this >= Q(Other);
	}

	bool operator>=(float Other) const
	{
		return *this >= Q(Other);
	}

	bool operator<(int Other) const
	{
		return *this < Q(Other);
	}

	bool operator<(float Other) const
	{
		return *this < Q(Other);
	}

	bool operator<=(int Other) const
	{
		return *this <= Q(Other);
	}

	bool operator<=(float Other) const
	{
		return *this <= Q(Other);
	}

	// Convenient arithmetic
	Q operator+(int Other) const
	{
		return *this + Q(Other);
	}

	Q operator+(float Other) const
	{
		return *this + Q(Other);
	}

	Q& operator+=(int Other)
	{
		return *this += Q(Other);
	}

	Q& operator+=(float Other)
	{
		return *this += Q(Other);
	}

	Q operator-(int Other) const
	{
		return *this - Q(Other);
	}

	Q operator-(float Other) const
	{
		return *this - Q(Other);
	}

	Q& operator-=(int Other)
	{
		return *this -= Q(Other);
	}

	Q& operator-=(float Other)
	{
		return *this -= Q(Other);
	}

	Q operator*(int Other) const
	{
		return *this * Q(Other);
	}

	Q operator*(float Other) const
	{
		return *this * Q(Other);
	}

	Q& operator*=(int Other)
	{
		return *this *= Q(Other);
	}

	Q& operator*=(float Other)
	{
		return *this *= Q(Other);
	}

	Q operator/(int Other) const
	{
		return *this / Q(Other);
	}

	Q operator/(float Other) const
	{
		return *this / Q(Other);
	}

	Q& operator/=(int Other)
	{
		return *this /= Q(Other);
	}

	Q& operator/=(float Other)
	{
		return *this /= Q(Other);
	}

	// Conversion / visualization
	float ToFloat() const
	{
		return ((float)Value) * Precision;
	}

	BaseType ToInt() const
	{
		return Value >> FractionalBits;
	}

	// General, algebraic and trigonometric functions
	static Q Negate(Q S)
	{
		return FromRaw(-S.Value);
	}
	
	static Q Abs(Q S)
	{
		if (S < Q(0.0f))
		{
			return Negate(S);
		}
		return S;
	}

	static Q Floor(Q S)
	{
		if (S < Q(0))
		{
			Q Pos = Negate(S);
			Q PosNoFrac = FromRaw(Pos.Value & IntegerMask);
			if (Pos == PosNoFrac)
			{
				return S;
			}
			return Negate(PosNoFrac + Q(1));
		}
		return FromRaw(S.Value & IntegerMask);
	}

	static Q Ceil(Q S)
	{
		if (S < Q(0))
		{
			Q Pos = Negate(S);
			Q PosNoFrac = FromRaw(Pos.Value & IntegerMask);
			return Negate(PosNoFrac);
		}
		Q NoFrac = FromRaw(S.Value & IntegerMask);
		if (NoFrac == S)
		{
			return S;
		}
		return NoFrac + Q(1);
	}

	static Q Truncate(Q S)
	{
		if (S < Q(0))
		{
			Q Pos = Negate(S);
			return Negate(FromRaw(Pos.Value & IntegerMask));
		}
		return FromRaw(S.Value & IntegerMask);
	}

	static Q ReciprocalSqrt(Q S)
	{
		return Q(1) / Sqrt(S);
	}

	static Q Sqrt(Q S)
	{
		if (S == Q(0))
		{
			return Q(0);
		}
		
		// Initial guess
		Q Res = FromRaw(S.Value >> 1); // S / 2

		// Five steps of the babylonian method
		for (int i = 0; i < 5; i++)
		{
			Res = FromRaw((Res + S / Res).Value >> 1); // Shift divides by 2
		}

		return Res;
	}

	static Q Fmod(Q X, Q Y)
	{
		return X - Floor(X / Y) * Y;
	}

	static Q Exp(Q S)
	{
		static const Q E = Q(M_E);

		bool Negate = false;
		if (S < Q(0))
		{
			Negate = true;
			S = Q::Negate(S);
		}

		Q Integral = IntPow(E, (unsigned int)S.ToInt());

		Q Frac;
		BaseType RawSFrac = S.Value & FractionalMask;
		if (RawSFrac > 0)
		{
			Frac = ExpTaylor(10, FromRaw(RawSFrac));
		}
		else
		{
			Frac = Q(1);
		}

		Q Total = Integral * Frac;

		return Negate ? Q(1) / Total : Total;
	}

private:
	// Taylor series sum
	// Taken from: https://www.geeksforgeeks.org/program-to-efficiently-calculate-ex/
	static Q ExpTaylor(unsigned int Iterations, Q S)
	{
		Q Sum = Q(1);

		for (unsigned int i = Iterations; i > 0; i--)
			Sum = Q(1) + S * Sum / Q((int)i);

		return Sum;
	}

	static Q CosApprox(Q Angle)
	{
		const Q C1 = Q(0.99940307f);
		const Q C2 = Q(-0.49558072f);
		const Q C3 = Q(0.03679168f);
		Q X2 = Angle * Angle;
		return (C1 + X2 * (C2 + C3 * X2));
	}

	// Taken from: http://www.ganssle.com/approx.htm
	// Note: this actually calculates tan(x * pi/4)
	// The main calling function adjusts for that.
	static Q TanApprox(Q Angle)
	{
		const Q C1 = Q(-3.6112171);
		const Q C2 = Q(-4.6133253);
		Q X2 = Angle * Angle;

		return (Angle * C1 / (C2 + X2));
	}

public:
	// Taken from: https://github.com/divideconcept/FastTrigo/blob/master/fasttrigo.cpp
	// Angle in radians
	static Q Cos(Q Angle)
	{
		static const Q Pi = Q(M_PI);
		static const Q HalfPi = Pi / Q(2);
		static const Q QuarterPi = HalfPi / Q(2);
		static const Q ThreeHalfPi = Q(3) * HalfPi;
		static const Q TwoPi = Q(2) * Pi;
		static const Q InvTwoPi = Q(1) / TwoPi;

		// Clamp angle to 0..2pi
		Angle = Angle - Floor(Angle * InvTwoPi) * TwoPi;
		Angle = Angle > Q(0) ? Angle : -Angle;

		if (Angle < HalfPi) return CosApprox(Angle);
		if (Angle < Pi) return -CosApprox(Pi - Angle);
		if (Angle < ThreeHalfPi) return -CosApprox(Angle - Pi);
		return CosApprox(TwoPi - Angle);
	}

	static Q Tan(Q Angle)
	{
		static const Q Pi = Q(M_PI);
		static const Q HalfPi = Pi / Q(2);
		static const Q QuarterPi = HalfPi / Q(2);
		static const Q ThreeHalfPi = Q(3) * HalfPi;
		static const Q TwoPi = Q(2) * Pi;
		static const Q InvTwoPi = Q(1) / TwoPi;
		static const Q FourOverPi = Q(4) / Pi;

		// Clamp angle to 0..2pi
		Angle = Angle - Floor(Angle * InvTwoPi) * TwoPi;
		Angle = Angle > Q(0) ? Angle : -Angle;

		BaseType Octant = (Angle / QuarterPi).ToInt();
		switch (Octant)
		{
		case 0: return TanApprox(Angle * FourOverPi);
		case 1: return Q(1) / TanApprox((HalfPi - Angle) * FourOverPi);
		case 2: return Q(-1) / TanApprox((Angle - HalfPi) * FourOverPi);
		case 3: return -TanApprox((Pi - Angle) * FourOverPi);
		case 4: return TanApprox((Angle - Pi) * FourOverPi);
		case 5: return Q(1) / TanApprox((ThreeHalfPi - Angle) * FourOverPi);
		case 6: return Q(-1) / TanApprox((Angle - ThreeHalfPi) * FourOverPi);
		case 7: return -TanApprox((TwoPi - Angle) * FourOverPi);
		}
		return FromRaw(0xDEADBEEF); // Cannot be reached, screw you compiler.
	}

	// Angle in radians
	static Q Sin(Q Angle)
	{
		static const Q Pi = Q(M_PI);
		static const Q HalfPi = Pi / Q(2);

		return Cos(HalfPi - Angle);
	}

	// Taken from: https://github.com/divideconcept/FastTrigo/blob/master/fasttrigo.cpp
	static Q Atan(Q Angle)
	{
		static const Q Pi = Q(M_PI);
		static const Q HalfPi = Pi / Q(2);
		static const Q QuarterPi = HalfPi / Q(2);

		return QuarterPi * Angle - Angle * (Abs(Angle) - Q(1)) * (Q(0.2447f) + Q(0.0663f) * Abs(Angle));
	}

	// Taken from: https://github.com/divideconcept/FastTrigo/blob/master/fasttrigo.cpp
	static Q Atan2(Q Y, Q X)
	{
		static const Q Pi = Q(M_PI);
		static const Q HalfPi = Pi / Q(2);

		if (Abs(X) > Abs(Y)) 
		{
			Q Atan = Q::Atan(Y / X);
			if (X > Q(0))
			{
				return Atan;
			}
			else
			{
				return Y > Q(0) ? Atan + Pi : Atan - Pi;
			}
		}
		else 
		{
			Q Atan = Q::Atan(X / Y);
			if (X > Q(0))
			{
				return Y > Q(0) ? HalfPi - Atan : -HalfPi - Atan;
			}
			else
			{
				return Y > Q(0) ? HalfPi + Atan : -HalfPi + Atan;
			}
		}
	}

	static Q Asin(Q X)
	{
		return Atan2(X, Sqrt((Q(1) + X) * (Q(1) - X)));
	}

	static Q Acos(Q X)
	{
		return Atan2(Sqrt((Q(1) + X) * (Q(1) - X)), X);
	}

	static Q IntPow(Q Base, unsigned int Exp)
	{
		Q Res = Q(1);

		for (;;)
		{
			if (Exp & 1)
				Res = Res * Base;
			Exp >>= 1;
			if (!Exp)
				break;
			Base = Base * Base;
		}

		return Res;
	}

	static Q Pow(Q Base, Q Exp)
	{
		bool Negate = false;
		if (Exp < Q(0))
		{
			Negate = true;
			Exp = Q::Negate(Exp);
		}

		// Divide our problem into integral and fractional parts.
		// Integral is easily solved.
		Q IntegralPow = IntPow(Base, (unsigned int)Exp.ToInt());

		// For the fractional part, we use pow(x, y) = exp(y * log(x))
		Q Fractional = FromRaw(Exp.Value & FractionalMask);
		Q FractionalPow = Q::Exp(Fractional * Log(Base));

		Q Total = IntegralPow * FractionalPow;

		return Negate ? Q(1) / Total : Total;
	}

	static Q Lerp(Q A, Q B, Q Alpha)
	{
		return A + Alpha * (B - A);
	}

public:
	// Idea taken from: https://stackoverflow.com/questions/4657468/fast-fixed-point-pow-log-exp-and-sqrt
	// We use the fact that Log(S) = Log(2) + Log(S/2)
	// And that Log(S) = Log(0.5) + Log(2 * S)
	// And also that Log(2) = -Log(0.5)
	// alongside a table with log values between 1 and 2 to get an acceptable approximation of our result.
	static Q Log(Q S)
	{
		static Q Log2 = Q(0.6931471805599453);
		static double const LogTable1to2[33] = {
			0.0,
			3.0771658666753687e-2,
			6.062462181643484e-2,
			8.961215868968714e-2,
			0.11778303565638346,
			0.1451820098444979,
			0.17185025692665923,
			0.19782574332991987,
			0.22314355131420976,
			0.24783616390458127,
			0.27193371548364176,
			0.2954642128938359,
			0.3184537311185346,
			0.3409265869705932,
			0.3629054936893685,
			0.38441169891033206,
			0.4054651081081644,
			0.4260843953109001,
			0.44628710262841953,
			0.46608972992459924,
			0.4855078157817008,
			0.5045560107523953,
			0.5232481437645479,
			0.5415972824327444,
			0.5596157879354227,
			0.5773153650348236,
			0.5947071077466928,
			0.6118015411059929,
			0.6286086594223741,
			0.6451379613735847,
			0.661398482245365,
			0.6773988235918061,
			0.6931471805599453
		};

		if (S <= Q(0))
		{
			// Lowest value possible.
			return FromRaw(0x8000000000000000LL);
		}

		int Shifts = 0;
		Q Normalized = S;
		if (S >= Q(2))
		{
			BaseType Integral = S.ToInt();
			while (Integral >= 2)
			{
				Integral >>= 1;
				Shifts++;
			}
			Normalized = FromRaw(S.Value >> Shifts);
		}
		else if (S < Q(1))
		{
			BaseType Fractional = S.Value & FractionalMask;
			while ((Fractional & IntegerMask) == 0)
			{
				Fractional <<= 1;
				Shifts--;
			}
			Normalized = FromRaw(Fractional);
		}

		BaseType Fractional = Normalized.Value & FractionalMask;
		int IntTableIndex = (int) (Fractional >> (FractionalBits - 5));
		Q Alpha = FromRaw((Fractional << 5) & FractionalMask);

		return Log2 * Q(Shifts) + Lerp(Q(LogTable1to2[IntTableIndex]), Q(LogTable1to2[IntTableIndex + 1]), Alpha);
	}
};

