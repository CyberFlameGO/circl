package group_test

import (
	"bytes"
	"crypto/rand"
	"fmt"
	"testing"

	"github.com/cloudflare/circl/group"
	"github.com/cloudflare/circl/internal/test"
)

var allGroups = []group.Group{
	group.P256,
	group.P384,
	group.P521,
	group.Ristretto255,
}

func TestGroup(t *testing.T) {
	const testTimes = 1 << 7
	for _, g := range allGroups {
		g := g
		n := g.(fmt.Stringer).String()
		t.Run(n+"/Add", func(tt *testing.T) { testAdd(tt, testTimes, g) })
		t.Run(n+"/Neg", func(tt *testing.T) { testNeg(tt, testTimes, g) })
		t.Run(n+"/Mul", func(tt *testing.T) { testMul(tt, testTimes, g) })
		t.Run(n+"/MulGen", func(tt *testing.T) { testMulGen(tt, testTimes, g) })
		t.Run(n+"/Order", func(tt *testing.T) { testOrder(tt, testTimes, g) })
		t.Run(n+"/Marshal", func(tt *testing.T) { testMarshal(tt, testTimes, g) })
		t.Run(n+"/Scalar", func(tt *testing.T) { testScalar(tt, testTimes, g) })
	}
}

func testAdd(t *testing.T, testTimes int, g group.Group) {
	Q := g.NewElement()
	for i := 0; i < testTimes; i++ {
		P := g.RandomElement(rand.Reader)

		got := Q.Dbl(P).Dbl(Q).Dbl(Q).Dbl(Q) // Q = 16P

		R := g.Identity()
		for j := 0; j < 16; j++ {
			R.Add(R, P)
		}
		want := R // R = 16P = P+P...+P
		if !got.IsEqual(want) {
			test.ReportError(t, got, want, P)
		}
	}
}

func testNeg(t *testing.T, testTimes int, g group.Group) {
	Q := g.NewElement()
	for i := 0; i < testTimes; i++ {
		P := g.RandomElement(rand.Reader)
		Q.Neg(P)
		Q.Add(Q, P)
		got := Q.IsIdentity()
		want := true
		if got != want {
			test.ReportError(t, got, want, P)
		}
	}
}

func testMul(t *testing.T, testTimes int, g group.Group) {
	Q := g.NewElement()
	kInv := g.NewScalar()
	for i := 0; i < testTimes; i++ {
		P := g.RandomElement(rand.Reader)
		k := g.RandomScalar(rand.Reader)
		kInv.Inv(k)

		Q.Mul(P, k)
		Q.Mul(Q, kInv)

		got := P
		want := Q
		if !got.IsEqual(want) {
			test.ReportError(t, got, want, P, k)
		}
	}
}

func testMulGen(t *testing.T, testTimes int, g group.Group) {
	G := g.Generator()
	P := g.NewElement()
	Q := g.NewElement()
	for i := 0; i < testTimes; i++ {
		k := g.RandomScalar(rand.Reader)

		P.Mul(G, k)
		Q.MulGen(k)

		got := P
		want := Q
		if !got.IsEqual(want) {
			test.ReportError(t, got, want, P, k)
		}
	}
}

func testOrder(t *testing.T, testTimes int, g group.Group) {
	Q := g.NewElement()
	order := g.Order()
	for i := 0; i < testTimes; i++ {
		P := g.RandomElement(rand.Reader)

		Q.Mul(P, order)
		got := Q.IsIdentity()
		want := true
		if got != want {
			test.ReportError(t, got, want, P)
		}
	}
}

func isZero(b []byte) bool {
	for i := 0; i < len(b); i++ {
		if b[i] != 0x00 {
			return false
		}
	}
	return true
}

func testMarshal(t *testing.T, testTimes int, g group.Group) {
	params := g.Params()
	I := g.Identity()
	got, err := I.MarshalBinary()
	test.CheckNoErr(t, err, "error on MarshalBinary")
	if !isZero(got) {
		test.ReportError(t, got, "Non-zero identity")
	}
	if l := uint(len(got)); !(l == 1 || l == params.ElementLength) {
		test.ReportError(t, l, params.ElementLength)
	}
	got, err = I.MarshalBinaryCompress()
	test.CheckNoErr(t, err, "error on MarshalBinaryCompress")
	if !isZero(got) {
		test.ReportError(t, got, "Non-zero identity")
	}
	if l := uint(len(got)); !(l == 1 || l == params.CompressedElementLength) {
		test.ReportError(t, l, params.CompressedElementLength)
	}
	II := g.NewElement()
	err = II.UnmarshalBinary(got)
	if err != nil || !I.IsEqual(II) {
		test.ReportError(t, I, II)
	}

	got1 := g.NewElement()
	got2 := g.NewElement()
	for i := 0; i < testTimes; i++ {
		x := g.RandomElement(rand.Reader)
		enc1, err1 := x.MarshalBinary()
		enc2, err2 := x.MarshalBinaryCompress()
		test.CheckNoErr(t, err1, "error on marshalling")
		test.CheckNoErr(t, err2, "error on marshalling compress")

		err1 = got1.UnmarshalBinary(enc1)
		err2 = got2.UnmarshalBinary(enc2)
		test.CheckNoErr(t, err1, "error on unmarshalling")
		test.CheckNoErr(t, err2, "error on unmarshalling compress")
		if !x.IsEqual(got1) {
			test.ReportError(t, got1, x)
		}
		if !x.IsEqual(got2) {
			test.ReportError(t, got2, x)
		}
		if l := uint(len(enc1)); l != params.ElementLength {
			test.ReportError(t, l, params.ElementLength)
		}
		if l := uint(len(enc2)); l != params.CompressedElementLength {
			test.ReportError(t, l, params.CompressedElementLength)
		}
	}
}

func testScalar(t *testing.T, testTimes int, g group.Group) {
	c := g.NewScalar()
	d := g.NewScalar()
	e := g.NewScalar()
	f := g.NewScalar()
	params := g.Params()
	for i := 0; i < testTimes; i++ {
		a := g.RandomScalar(rand.Reader)
		b := g.RandomScalar(rand.Reader)
		c.Add(a, b)
		d.Sub(a, b)
		e.Mul(c, d)

		c.Mul(a, a)
		d.Mul(b, b)
		d.Neg(d)
		f.Add(c, d)
		enc1, err1 := e.MarshalBinary()
		enc2, err2 := f.MarshalBinary()
		if err1 != nil || err2 != nil || !bytes.Equal(enc1, enc2) {
			test.ReportError(t, enc1, enc2, a, b)
		}
		if l := uint(len(enc1)); l != params.ScalarLength {
			test.ReportError(t, l, params.ScalarLength)
		}
	}
}

func BenchmarkElement(b *testing.B) {
	for _, g := range allGroups {
		x := g.RandomElement(rand.Reader)
		y := g.RandomElement(rand.Reader)
		n := g.RandomScalar(rand.Reader)
		name := g.(fmt.Stringer).String()
		b.Run(name+"/Add", func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				x.Add(x, y)
			}
		})
		b.Run(name+"/Dbl", func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				x.Dbl(x)
			}
		})
		b.Run(name+"/Mul", func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				y.Mul(x, n)
			}
		})
		b.Run(name+"/MulGen", func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				x.MulGen(n)
			}
		})
	}
}

func BenchmarkScalar(b *testing.B) {
	for _, g := range allGroups {
		x := g.RandomScalar(rand.Reader)
		y := g.RandomScalar(rand.Reader)
		name := g.(fmt.Stringer).String()
		b.Run(name+"/Add", func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				x.Add(x, y)
			}
		})
		b.Run(name+"/Mul", func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				x.Mul(x, y)
			}
		})
		b.Run(name+"/Inv", func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				y.Inv(x)
			}
		})
	}
}
