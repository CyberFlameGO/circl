package ff

import "fmt"

// Fp12Alt represents elements of Fp4[w]/w^3-t
type Fp12Alt [3]Fp4

// LineValue a represents a[0]+a[1]*w^2+a[2]*w^3, with all values in Fp2.
// This lets us shave off a number of Fp2 multiplications.
type LineValue [3]Fp2

func (z Fp12Alt) String() string { return fmt.Sprintf("%s + ( %s )*w + ( %s )*w^2", z[0], z[1], z[2]) }
func (z Fp12Alt) IsEqual(x *Fp12Alt) int {
	return z[0].IsEqual(&x[0]) & z[1].IsEqual(&x[1]) & z[2].IsEqual(&x[2])
}
func (z *Fp12Alt) FromFp12(x *Fp12) {
	// To understand this function, write everything as w^6-(u+1).
	// v = w^2
	// t = w^3
	z[0][0] = x[0][0] // w^0
	z[1][0] = x[1][0] // w^1
	z[2][0] = x[0][1] // w^2
	z[0][1] = x[1][1] // w^3
	z[1][1] = x[0][2] // w^4
	z[2][1] = x[1][2] // w^5
}

func (z *Fp12) FromFp12Alt(x *Fp12Alt) {
	z[0][0] = x[0][0] // w^0
	z[1][0] = x[1][0] // w^1
	z[0][1] = x[2][0] // w^2
	z[1][1] = x[0][1] // w^3
	z[0][2] = x[1][1] // w^4
	z[1][2] = x[2][1] // w^5
}

func (z *Fp12Alt) Add(x *Fp12Alt, y *Fp12Alt) {
	z[0].Add(&x[0], &y[0])
	z[1].Add(&x[1], &y[1])
	z[2].Add(&x[2], &y[2])
}

func (z *Fp12Alt) SetOne() {
	z[0].SetOne()
	z[1] = Fp4{}
	z[2] = Fp4{}
}

func (z *Fp12Alt) Mul(x *Fp12Alt, y *Fp12Alt) {
	// See Multiplication and Squaring in Pairing Friendly Fields
	var v0, v1, v2, p0, p1, p2, tx, ty Fp4
	v0.Mul(&x[0], &y[0])
	v1.Mul(&x[1], &y[1])
	v2.Mul(&x[2], &y[2])

	tx.Add(&x[1], &x[2])
	ty.Add(&y[1], &y[2])
	p0.Mul(&tx, &ty)

	tx.Add(&x[0], &x[1])
	ty.Add(&y[0], &y[1])
	p1.Mul(&tx, &ty)

	tx.Add(&x[0], &x[2])
	ty.Add(&y[0], &y[2])
	p2.Mul(&tx, &ty)

	z[0].Sub(&p0, &v1)
	z[0].Sub(&z[0], &v2)
	z[0].MulT(&z[0])
	z[0].Add(&z[0], &v0)

	z[1].Sub(&p1, &v0)
	z[1].Sub(&z[1], &v1)
	tx.MulT(&v2)
	z[1].Add(&z[1], &tx)

	z[2].Sub(&p2, &v0)
	z[2].Add(&z[2], &v1)
	z[2].Sub(&z[2], &v2)
}

func (z *Fp12Alt) Sqr(x *Fp12Alt) {
	var s0, s1, s2, s3, s4, t Fp4
	s0.Sqr(&x[0])
	s1.Mul(&x[0], &x[1])
	s1.Add(&s1, &s1)
	t.Add(&x[0], &x[2])
	t.Sub(&t, &x[1])
	s2.Sqr(&t)
	s3.Mul(&x[1], &x[2])
	s3.Add(&s3, &s3)
	s4.Sqr(&x[2])

	z[0].MulT(&s3)
	z[0].Add(&z[0], &s0)

	z[1].MulT(&s4)
	z[1].Add(&z[1], &s1)

	z[2].Add(&s1, &s2)
	z[2].Add(&z[2], &s3)
	z[2].Sub(&z[2], &s0)
	z[2].Sub(&z[2], &s4)
}

func (z *Fp12Alt) MulLine(x *Fp12Alt, y *LineValue) {
	var y0 Fp4
	var y2 Fp2
	var v0, v2, p0, p1, p2, tx, ty Fp4

	y0[0] = y[0]
	y0[1] = y[2]
	y2 = y[1]

	v0.Mul(&x[0], &y0)
	v2.MulSubfield(&x[2], &y2)

	tx.Add(&x[1], &x[2])
	p0.MulSubfield(&tx, &y2)

	tx.Add(&x[0], &x[1])
	p1.Mul(&tx, &y0)

	tx.Add(&x[0], &x[2])
	ty = y0
	ty[0].Add(&ty[0], &y2)
	p2.Mul(&tx, &ty)

	z[0].Sub(&p0, &v2)
	z[0].MulT(&z[0])
	z[0].Add(&z[0], &v0)

	z[1].Sub(&p1, &v0)
	tx.MulT(&v2)
	z[1].Add(&z[1], &tx)

	z[2].Sub(&p2, &v0)
	z[2].Sub(&z[2], &v2)
}

func (z *LineValue) IsZero() int {
	return z[0].IsZero() & z[1].IsZero() & z[2].IsZero()
}

func (z *LineValue) SetOne() {
	z[0].SetOne()
	z[1] = Fp2{}
	z[2] = Fp2{}
}
